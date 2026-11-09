"""
qPCR-Auswertung für Bio-Rad CFX (384-Well)
==========================================
Start:  streamlit run app.py
"""
import io
import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

try:
    from scipy import stats as _stats
except Exception:  # scipy optional
    _stats = None

# --------------------------------------------------------------------------------------
# Konstanten
# --------------------------------------------------------------------------------------
ROWS = list("ABCDEFGHIJKLMNOP")
NCOLS = 24
ALL_WELLS = [f"{r}{c}" for r in ROWS for c in range(1, NCOLS + 1)]
LABEL_COLS = ["Gen", "Treatment", "Probe", "Replikat", "Ausschluss"]

# Kategorische Farben (feste Reihenfolge) + Sequentielle Blau-Rampe
CAT_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
              "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ_BLUES = [[0.0, "#0d366b"], [0.25, "#1c5cab"], [0.5, "#3987e5"],
             [0.75, "#86b6ef"], [1.0, "#cde2fb"]]  # dunkel = niedriger Cq
EMPTY_COLOR = "#e4e3df"
NA_COLOR = "#b9b8b2"

st.set_page_config(page_title="qPCR-Auswertung · Bio-Rad 384", layout="wide",
                   page_icon="🧬")


# --------------------------------------------------------------------------------------
# Datei einlesen
# --------------------------------------------------------------------------------------
def norm_well(w):
    m = re.match(r"^\s*([A-Pa-p])\s*0*(\d{1,2})\s*$", str(w))
    if not m:
        return None
    c = int(m.group(2))
    if not 1 <= c <= NCOLS:
        return None
    return f"{m.group(1).upper()}{c}"


def _is_header(cells):
    cells = [str(x).strip().lower() for x in cells]
    return "well" in cells and any(c in ("cq", "ct", "c(t)", "cт") for c in cells)


def _to_num(s):
    return pd.to_numeric(s.astype(str).str.strip().str.replace(",", ".", regex=False),
                         errors="coerce")


@st.cache_data(show_spinner=False)
def load_file(name: str, data: bytes) -> pd.DataFrame:
    df = None
    if name.lower().endswith((".xlsx", ".xlsm", ".xls")):
        sheets = pd.read_excel(io.BytesIO(data), header=None, dtype=str, sheet_name=None)
        for _, raw in sheets.items():
            raw = raw.fillna("")
            for i in range(min(len(raw), 60)):
                if _is_header(raw.iloc[i].tolist()):
                    df = raw.iloc[i + 1:].copy()
                    df.columns = [str(c).strip() for c in raw.iloc[i].tolist()]
                    break
            if df is not None:
                break
    else:
        text = None
        for enc in ("utf-8-sig", "cp1252", "latin-1"):
            try:
                text = data.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        lines = text.splitlines()
        hdr = None
        for i, line in enumerate(lines[:80]):
            for sep in (";", "\t", ","):
                if _is_header(line.split(sep)):
                    hdr, use_sep = i, sep
                    break
            if hdr is not None:
                break
        if hdr is not None:
            df = pd.read_csv(io.StringIO("\n".join(lines[hdr:])), sep=use_sep, dtype=str,
                             keep_default_na=False, engine="python")
            df.columns = [str(c).strip().strip('"') for c in df.columns]

    if df is None:
        raise ValueError("Keine Kopfzeile mit 'Well' und 'Cq' gefunden. "
                         "Ist das ein Bio-Rad CFX Export (Quantification Cq Results)?")

    lower = {c.lower(): c for c in df.columns if c}
    pick = lambda *names: next((lower[n] for n in names if n in lower), None)  # noqa: E731
    out = pd.DataFrame({"Well": df[pick("well")].map(norm_well)})
    out["Cq"] = _to_num(df[pick("cq", "ct", "c(t)")])
    for target, names in {"Fluor": ("fluor", "dye"), "Target": ("target",),
                          "Sample": ("sample",), "Content": ("content",)}.items():
        c = pick(*names)
        out[target] = df[c].astype(str).str.strip().replace({"nan": ""}) if c else ""
    out = out[out["Well"].notna()].reset_index(drop=True)
    if out.empty:
        raise ValueError("Keine gültigen Well-Bezeichnungen (A1–P24) gefunden.")
    return out


# --------------------------------------------------------------------------------------
# Hilfsfunktionen Wells / Bereiche
# --------------------------------------------------------------------------------------
def rc(well):
    return ROWS.index(well[0]), int(well[1:])


def sort_wells(wells):
    return sorted(wells, key=lambda w: rc(w))


def _parse_ref(ref):
    m = re.match(r"^([A-Pa-p])?(\d{1,2})?$", ref)
    if not m or (m.group(1) is None and m.group(2) is None):
        raise ValueError(f"Ungültige Angabe: '{ref}'")
    r = ROWS.index(m.group(1).upper()) if m.group(1) else None
    c = int(m.group(2)) if m.group(2) else None
    if c is not None and not 1 <= c <= NCOLS:
        raise ValueError(f"Spalte außerhalb 1–24: '{ref}'")
    return r, c


def parse_range(text):
    """'A1:C16; E1-E16, B, 5:7' -> Liste von Wells."""
    text = re.sub(r"\s*([:\-])\s*", r"\1", text.strip())
    wells = set()
    for part in re.split(r"[;,\s]+", text):
        if not part:
            continue
        a, b = (re.split(r"[:\-]", part, maxsplit=1) + [None])[:2]
        b = b or a
        (ra, ca), (rb, cb) = _parse_ref(a), _parse_ref(b)
        rr = [x for x in (ra, rb) if x is not None]
        cc = [x for x in (ca, cb) if x is not None]
        r0, r1 = (min(rr), max(rr)) if rr else (0, len(ROWS) - 1)
        c0, c1 = (min(cc), max(cc)) if cc else (1, NCOLS)
        for r in range(r0, r1 + 1):
            for c in range(c0, c1 + 1):
                wells.add(f"{ROWS[r]}{c}")
    return sort_wells(wells)


def describe_wells(wells):
    if not wells:
        return "keine"
    rows = sorted({rc(w)[0] for w in wells})
    cols = sorted({rc(w)[1] for w in wells})
    box = (len(rows) * len(cols) == len(wells)
           and rows == list(range(rows[0], rows[-1] + 1))
           and cols == list(range(cols[0], cols[-1] + 1)))
    span = f"{ROWS[rows[0]]}{cols[0]}:{ROWS[rows[-1]]}{cols[-1]}"
    return f"{len(wells)} Wells ({span}{'' if box else ', nicht rechteckig'})"


def sample_blocks(wells, n_rep, rep_dir, order, start):
    """Teilt Wells in Proben-Blöcke à n_rep Replikate.
    rep_dir: 'h' = Replikate nebeneinander in einer Zeile, 'v' = untereinander in einer Spalte
    order:   'down' = Probennummern laufen zuerst nach unten, 'right' = zuerst nach rechts
    """
    coords = [(*rc(w), w) for w in wells]
    chunks = []
    if rep_dir == "h":
        for r in sorted({c[0] for c in coords}):
            line = sorted([c for c in coords if c[0] == r], key=lambda x: x[1])
            for i in range(0, len(line), n_rep):
                ch = line[i:i + n_rep]
                chunks.append(((r, ch[0][1]), [x[2] for x in ch]))
    else:
        for col in sorted({c[1] for c in coords}):
            line = sorted([c for c in coords if c[1] == col], key=lambda x: x[0])
            for i in range(0, len(line), n_rep):
                ch = line[i:i + n_rep]
                chunks.append(((ch[0][0], col), [x[2] for x in ch]))
    if order == "down":
        chunks.sort(key=lambda x: (x[0][1], x[0][0]))
    else:
        chunks.sort(key=lambda x: (x[0][0], x[0][1]))
    result, incomplete = {}, 0
    for k, (_, ws) in enumerate(chunks):
        if len(ws) < n_rep:
            incomplete += 1
        for j, w in enumerate(ws):
            result[w] = (start + k, j + 1)
    return result, len(chunks), incomplete


# --------------------------------------------------------------------------------------
# Session-State
# --------------------------------------------------------------------------------------
def empty_labels(wells):
    return pd.DataFrame({"Gen": "", "Treatment": "", "Probe": pd.Series([pd.NA] * len(wells), dtype="Int64"),
                         "Replikat": pd.Series([pd.NA] * len(wells), dtype="Int64"),
                         "Ausschluss": False}, index=pd.Index(wells, name="Well"))


def clean_labels(lab):
    lab = lab.copy()
    for c in ("Gen", "Treatment"):
        lab[c] = lab[c].fillna("").astype(str).str.strip().replace({"<NA>": "", "nan": "", "None": ""})
    for c in ("Probe", "Replikat"):
        lab[c] = pd.to_numeric(lab[c], errors="coerce").round().astype("Int64")
    lab["Ausschluss"] = lab["Ausschluss"].fillna(False).astype(bool)
    return lab[LABEL_COLS]


def push_history():
    hist = st.session_state.setdefault("history", [])
    hist.append(st.session_state.labels.copy())
    del hist[:-30]


def set_labels(new):
    push_history()
    st.session_state.labels = clean_labels(new)
    st.session_state.ver += 1  # Editor zurücksetzen


def init_state(raw, sig):
    if st.session_state.get("sig") == sig:
        return
    wells = sort_wells(raw["Well"].unique().tolist())
    lab = empty_labels(wells)
    # Vorbelegung aus der Datei, falls dort schon Target/Sample stehen
    first = raw.drop_duplicates("Well").set_index("Well")
    lab["Gen"] = first["Target"].reindex(wells).fillna("").astype(str)
    st.session_state.update(sig=sig, labels=clean_labels(lab), history=[], ver=0,
                            sel=[], last_evt=None, plot_ver=0)


# --------------------------------------------------------------------------------------
# Plattenansicht
# --------------------------------------------------------------------------------------
def plate_figure(d, color_by, show_text, selected, height=560):
    """d: DataFrame (Index Well) mit Cq + Labels."""
    fig = go.Figure()
    # Leere Positionen (nicht in der Datei)
    absent = [w for w in ALL_WELLS if w not in d.index]
    if absent:
        fig.add_trace(go.Scatter(
            x=[rc(w)[1] for w in absent], y=[rc(w)[0] for w in absent], mode="markers",
            marker=dict(symbol="circle-open", size=14, color=EMPTY_COLOR, line=dict(width=1)),
            hoverinfo="skip", showlegend=False, name="leer"))

    def cdata(sub):
        return np.stack([sub.index, sub["Cq"].map(lambda v: "—" if pd.isna(v) else f"{v:.2f}"),
                         sub["Gen"], sub["Treatment"],
                         sub["Probe"].astype(str).replace("<NA>", ""),
                         sub["Replikat"].astype(str).replace("<NA>", ""),
                         np.where(sub["Ausschluss"], "ja", "nein")], axis=-1)

    hover = ("<b>%{customdata[0]}</b>  Cq %{customdata[1]}<br>Gen: %{customdata[2]}"
             "<br>Treatment: %{customdata[3]}<br>Probe: %{customdata[4]}  ·  Rep: %{customdata[5]}"
             "<br>Ausgeschlossen: %{customdata[6]}<extra></extra>")
    base_marker = dict(symbol="square", size=21, line=dict(width=1.5, color="#fcfcfb"))

    def text_for(sub):
        if not show_text:
            return None
        if color_by == "Cq":
            return sub["Cq"].map(lambda v: "" if pd.isna(v) else f"{v:.1f}")
        if color_by == "Probe":
            return sub["Probe"].astype(str).replace("<NA>", "")
        if color_by == "Replikat":
            return sub["Replikat"].astype(str).replace("<NA>", "")
        return sub[color_by].str[:3]

    if color_by in ("Cq", "Probe"):
        col = "Cq" if color_by == "Cq" else "Probe"
        has = d[d[col].notna()]
        nan = d[d[col].isna()]
        if len(has):
            vals = has[col].astype(float)
            fig.add_trace(go.Scatter(
                x=[rc(w)[1] for w in has.index], y=[rc(w)[0] for w in has.index],
                mode="markers+text" if show_text else "markers", text=text_for(has),
                textfont=dict(size=8, color=list(np.where(
                    (vals - vals.min()) / max(vals.max() - vals.min(), 1e-9) < 0.5,
                    "#ffffff" if col == "Cq" else "#0b0b0b",
                    "#0b0b0b" if col == "Cq" else "#ffffff"))),
                marker=dict(**base_marker, color=vals,
                            colorscale=SEQ_BLUES if col == "Cq" else SEQ_BLUES[::-1],
                            colorbar=dict(title=col, thickness=12, len=0.8)),
                customdata=cdata(has), hovertemplate=hover, showlegend=False, name=col))
        if len(nan):
            fig.add_trace(go.Scatter(
                x=[rc(w)[1] for w in nan.index], y=[rc(w)[0] for w in nan.index], mode="markers",
                marker=dict(**{**base_marker, "color": NA_COLOR}),
                customdata=cdata(nan), hovertemplate=hover,
                name="kein Cq" if col == "Cq" else "ohne Probe"))
    else:
        cats = [c for c in pd.unique(d[color_by].astype(str)) if c not in ("", "<NA>")]
        cats = sorted(cats, key=lambda x: (len(x), x)) if color_by == "Replikat" else cats
        for i, cat in enumerate(cats):
            sub = d[d[color_by].astype(str) == cat]
            # >8 Kategorien: Farben wiederholen, Identität über Text/Hover (sekundäre Kodierung)
            color = CAT_COLORS[i % len(CAT_COLORS)]
            fig.add_trace(go.Scatter(
                x=[rc(w)[1] for w in sub.index], y=[rc(w)[0] for w in sub.index],
                mode="markers+text" if show_text else "markers", text=text_for(sub),
                textfont=dict(size=8, color="#0b0b0b"),
                marker=dict(**{**base_marker, "color": color}),
                customdata=cdata(sub), hovertemplate=hover, name=cat))
        rest = d[~d[color_by].astype(str).isin(cats)]
        if len(rest):
            fig.add_trace(go.Scatter(
                x=[rc(w)[1] for w in rest.index], y=[rc(w)[0] for w in rest.index], mode="markers",
                marker=dict(**{**base_marker, "color": NA_COLOR}),
                customdata=cdata(rest), hovertemplate=hover, name="unbeschriftet"))

    exc = d[d["Ausschluss"]]
    if len(exc):
        fig.add_trace(go.Scatter(
            x=[rc(w)[1] for w in exc.index], y=[rc(w)[0] for w in exc.index], mode="markers",
            marker=dict(symbol="x-thin-open", size=16, line=dict(width=2.5, color="#0b0b0b")),
            hoverinfo="skip", name="ausgeschlossen"))
    if selected:
        fig.add_trace(go.Scatter(
            x=[rc(w)[1] for w in selected], y=[rc(w)[0] for w in selected], mode="markers",
            marker=dict(symbol="square-open", size=24, line=dict(width=2.5, color="#e34948")),
            hoverinfo="skip", name="Auswahl"))

    fig.update_layout(
        height=height, margin=dict(l=10, r=10, t=30, b=10), dragmode="select",
        plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)", clickmode="event+select",
        legend=dict(orientation="h", y=-0.03, yanchor="top", x=0),
        xaxis=dict(range=[0.4, NCOLS + 0.6], tickvals=list(range(1, NCOLS + 1)), side="top",
                   showgrid=False, zeroline=False, fixedrange=True),
        yaxis=dict(range=[len(ROWS) - 0.4, -0.6], tickvals=list(range(len(ROWS))), ticktext=ROWS,
                   showgrid=False, zeroline=False, fixedrange=True),
    )
    return fig


def plate_with_selection(d, color_by, show_text, key_prefix, height=560):
    fig = plate_figure(d, color_by, show_text, st.session_state.sel, height)
    evt = st.plotly_chart(fig, key=f"{key_prefix}_{st.session_state.plot_ver}", on_select="rerun",
                          selection_mode=("points", "box", "lasso"), use_container_width=True,
                          config={"displaylogo": False, "modeBarButtonsToRemove": ["zoom", "pan"]})
    pts = (evt or {}).get("selection", {}).get("points", []) if evt else []
    wells = set()
    for p in pts:
        try:
            w = f"{ROWS[int(round(p['y']))]}{int(round(p['x']))}"
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        if w in d.index:
            wells.add(w)
    wells = sort_wells(wells)
    if wells and wells != st.session_state.last_evt:
        st.session_state.last_evt = wells
        if st.session_state.get("sel_mode", "ersetzen") == "hinzufügen":
            st.session_state.sel = sort_wells(set(st.session_state.sel) | set(wells))
        else:
            st.session_state.sel = wells
        st.rerun()


def clear_selection():
    st.session_state.sel = []
    st.session_state.last_evt = None
    st.session_state.plot_ver += 1


# --------------------------------------------------------------------------------------
# Auswertung
# --------------------------------------------------------------------------------------
def replicate_summary(d, sd_thr):
    v = d[(~d["Ausschluss"]) & (d["Gen"] != "") & (d["Treatment"] != "") & d["Probe"].notna()]
    if v.empty:
        return pd.DataFrame()
    g = v.reset_index().groupby(["Gen", "Treatment", "Probe"], sort=True)
    s = g.agg(Wells=("Well", lambda x: ", ".join(sort_wells(x))), n_Wells=("Well", "size"),
              n_Cq=("Cq_eff", "count"), Cq_Mittel=("Cq_eff", "mean"), Cq_SD=("Cq_eff", "std"),
              Cq_Min=("Cq_eff", "min"), Cq_Max=("Cq_eff", "max")).reset_index()
    s["Flag"] = np.where(s["n_Cq"] == 0, "kein Cq",
                np.where(s["n_Cq"] < s["n_Wells"], "Replikat ohne Cq",
                np.where(s["Cq_SD"] > sd_thr, f"SD > {sd_thr}", "")))
    return s


def outlier_suggestions(d, sd_thr):
    """Für Gruppen mit ≥3 Cq-Werten und SD > Schwelle: Well mit größter Abweichung vom Median."""
    v = d[(~d["Ausschluss"]) & (d["Gen"] != "") & (d["Treatment"] != "")
          & d["Probe"].notna() & d["Cq_eff"].notna()].reset_index()
    sug = []
    for _, grp in v.groupby(["Gen", "Treatment", "Probe"]):
        if len(grp) >= 3 and grp["Cq_eff"].std() > sd_thr:
            dev = (grp["Cq_eff"] - grp["Cq_eff"].median()).abs()
            w = grp.loc[dev.idxmax(), "Well"]
            rest = grp[grp["Well"] != w]["Cq_eff"]
            if rest.std() <= sd_thr:  # nur vorschlagen, wenn es danach passt
                sug.append(w)
    return sort_wells(sug)


def ddcq(summary, refs, control):
    s = summary.dropna(subset=["Cq_Mittel"])
    ref = (s[s["Gen"].isin(refs)].groupby(["Treatment", "Probe"])
           .agg(Ref_Cq=("Cq_Mittel", "mean"), n_ref=("Gen", "nunique")).reset_index())
    ref = ref[ref["n_ref"] == len(refs)].drop(columns="n_ref")
    t = s[~s["Gen"].isin(refs)].merge(ref, on=["Treatment", "Probe"], how="inner")
    if t.empty:
        return pd.DataFrame(), pd.DataFrame()
    t["dCq"] = t["Cq_Mittel"] - t["Ref_Cq"]
    ctrl = t[t["Treatment"] == control].groupby("Gen")["dCq"].mean().rename("dCq_Kontrolle")
    t = t.merge(ctrl, on="Gen", how="left")
    t["ddCq"] = t["dCq"] - t["dCq_Kontrolle"]
    t["Fold_Change"] = 2.0 ** (-t["ddCq"])
    per = t[["Gen", "Treatment", "Probe", "Cq_Mittel", "Ref_Cq", "dCq", "ddCq", "Fold_Change"]]

    rows = []
    for (gen, tr), grp in per.groupby(["Gen", "Treatment"], sort=False):
        r = dict(Gen=gen, Treatment=tr, n_Proben=len(grp), dCq_Mittel=grp["dCq"].mean(),
                 ddCq_Mittel=grp["ddCq"].mean(), FC_Mittel=grp["Fold_Change"].mean(),
                 FC_SD=grp["Fold_Change"].std(), FC_SEM=grp["Fold_Change"].sem(),
                 FC_aus_ddCq=2.0 ** (-grp["ddCq"].mean()), p_Wert=np.nan)
        if _stats is not None and tr != control:
            c = per[(per["Gen"] == gen) & (per["Treatment"] == control)]["dCq"]
            if len(grp) >= 2 and len(c) >= 2:
                r["p_Wert"] = _stats.ttest_ind(grp["dCq"], c, equal_var=False).pvalue
        rows.append(r)
    return per, pd.DataFrame(rows)


def fc_figure(per, summ, control, log2):
    genes = list(dict.fromkeys(summ["Gen"]))
    treats = list(dict.fromkeys([control] + [t for t in summ["Treatment"] if t != control]))
    fig = go.Figure()
    for i, tr in enumerate(treats):
        s = summ[summ["Treatment"] == tr].set_index("Gen").reindex(genes)
        y = np.log2(s["FC_Mittel"]) if log2 else s["FC_Mittel"]
        color = CAT_COLORS[i % len(CAT_COLORS)]
        err = None if log2 else dict(type="data", array=s["FC_SD"], color="#52514e", thickness=1.2, width=4)
        fig.add_trace(go.Bar(x=genes, y=y, name=tr, marker=dict(color=color, line=dict(width=2, color="#fcfcfb")),
                             error_y=err, offsetgroup=tr,
                             hovertemplate=f"<b>%{{x}}</b> · {tr}<br>FC: %{{customdata[0]:.3f}}"
                                           "<br>n = %{customdata[1]}<extra></extra>",
                             customdata=np.stack([s["FC_Mittel"], s["n_Proben"]], axis=-1)))
        p = per[per["Treatment"] == tr]
        fig.add_trace(go.Box(x=p["Gen"], y=np.log2(p["Fold_Change"]) if log2 else p["Fold_Change"],
                             name=tr, offsetgroup=tr, boxpoints="all", jitter=0.4, pointpos=0,
                             fillcolor="rgba(0,0,0,0)", line=dict(color="rgba(0,0,0,0)"),
                             marker=dict(color="#0b0b0b", size=8, opacity=0.65,
                                         line=dict(width=1, color="#fcfcfb")),
                             showlegend=False, customdata=p["Probe"],
                             hovertemplate="Probe %{customdata}: %{y:.3f}<extra></extra>"))
    fig.update_layout(barmode="group", boxmode="group", height=460,
                      yaxis_title="log2 Fold Change" if log2 else "Fold Change (2^-ΔΔCq)",
                      margin=dict(l=10, r=10, t=30, b=10), legend_title="Treatment",
                      plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)")
    fig.update_yaxes(gridcolor="rgba(128,128,128,0.15)", zeroline=True,
                     zerolinecolor="rgba(128,128,128,0.5)")
    if not log2:
        fig.add_hline(y=1, line=dict(color="#52514e", width=1, dash="dot"))
    return fig


def to_excel(sheets: dict) -> bytes:
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        for name, df in sheets.items():
            if df is not None and len(df):
                df.to_excel(xw, sheet_name=name[:31], index=False)
    return buf.getvalue()


# ======================================================================================
# UI
# ======================================================================================
st.title("🧬 qPCR-Auswertung · Bio-Rad CFX 384")

with st.sidebar:
    st.header("1 · Daten")
    up = st.file_uploader("Bio-Rad Export (Quantification Cq Results)", type=["xlsx", "xls", "csv", "txt"])
    st.header("Einstellungen")
    cutoff_on = st.checkbox("Cq ≥ Grenzwert als 'nicht detektiert' werten", value=True)
    cutoff = st.number_input("Grenzwert Cq", value=40.0, step=0.5, disabled=not cutoff_on)
    sd_thr = st.number_input("Max. SD technischer Replikate", value=0.5, step=0.05, min_value=0.0)

if up is None:
    st.info("⬅️ Lade links eine Bio-Rad CFX-Exportdatei (.xlsx oder .csv) hoch. "
            "Erwartet wird die Tabelle **Quantification Cq Results** mit den Spalten *Well* und *Cq*.")
    st.stop()

try:
    raw = load_file(up.name, up.getvalue())
except Exception as e:  # noqa: BLE001
    st.error(f"Datei konnte nicht gelesen werden: {e}")
    st.stop()

fluors = [f for f in raw["Fluor"].unique() if f]
if len(fluors) > 1:
    fl = st.sidebar.selectbox("Fluorophor / Kanal", fluors)
    raw = raw[raw["Fluor"] == fl]
raw = raw.drop_duplicates("Well")
init_state(raw, (up.name, up.size, tuple(fluors)))

# Layout-Vorlage laden / speichern
with st.sidebar:
    st.header("Layout-Vorlage")
    st.caption("Beschriftung als CSV speichern und für die nächste Platte wiederverwenden.")
    st.download_button("💾 Layout speichern (CSV)",
                       st.session_state.labels.reset_index().to_csv(index=False, sep=";").encode("utf-8-sig"),
                       file_name="plattenlayout.csv", mime="text/csv", use_container_width=True)
    tpl = st.file_uploader("Layout laden", type=["csv"], key="tpl")
    if tpl is not None and st.button("Layout anwenden", use_container_width=True):
        try:
            t = pd.read_csv(tpl, sep=None, engine="python", dtype=str, encoding="utf-8-sig")
            t["Well"] = t["Well"].map(norm_well)
            t = t.dropna(subset=["Well"]).set_index("Well")
            new = st.session_state.labels.copy()
            common = new.index.intersection(t.index)
            for c in LABEL_COLS:
                if c in t.columns:
                    vals = t.loc[common, c]
                    if c == "Ausschluss":
                        vals = vals.astype(str).str.lower().isin(["true", "1", "ja", "x", "wahr"])
                    new.loc[common, c] = vals
            set_labels(new)
            st.success(f"Layout auf {len(common)} Wells angewendet.")
        except Exception as e:  # noqa: BLE001
            st.error(f"Layout konnte nicht gelesen werden: {e}")

# Arbeitsdaten
labels = st.session_state.labels
d = raw.set_index("Well")[["Cq"]].join(labels, how="left")
d["Cq_eff"] = d["Cq"].where(~(d["Cq"] >= cutoff)) if cutoff_on else d["Cq"]

n_lab = int(((labels["Gen"] != "") & (labels["Treatment"] != "") & labels["Probe"].notna()).sum())
c1, c2, c3, c4 = st.columns(4)
c1.metric("Wells in Datei", len(d))
c2.metric("mit Cq", int(d["Cq_eff"].notna().sum()))
c3.metric("vollständig beschriftet", n_lab)
c4.metric("ausgeschlossen", int(labels["Ausschluss"].sum()))

tab_plate, tab_label, tab_res = st.tabs(["📊 Platte", "🏷️ Beschriften", "📈 Auswertung"])

# --------------------------------------------------------------------------------------
with tab_plate:
    a, b = st.columns([3, 1])
    color_by = a.radio("Färben nach", ["Cq", "Gen", "Treatment", "Probe", "Replikat"], horizontal=True,
                       key="cb_plate")
    show_text = b.toggle("Werte in Wells anzeigen", value=True, key="txt_plate")
    view = d.copy()
    view["Cq"] = view["Cq_eff"]
    st.plotly_chart(plate_figure(view, color_by, show_text, []), use_container_width=True,
                    config={"displaylogo": False})
    st.caption("Grau = kein Cq (bzw. ≥ Grenzwert) · leere Kreise = Well nicht in der Datei · ✕ = ausgeschlossen. "
               "Dunkleres Blau = niedrigerer Cq.")
    with st.expander("Cq-Matrix (16 × 24)"):
        mat = pd.DataFrame(index=ROWS, columns=range(1, NCOLS + 1), dtype=float)
        for w, v in d["Cq_eff"].items():
            r, c = rc(w)
            mat.iloc[r, c - 1] = v
        mat = mat.astype(float)
        lo, hi = np.nanmin(mat.values) if mat.notna().any().any() else 0, np.nanmax(mat.values) if mat.notna().any().any() else 1

        def _cell_color(v):  # eigene Farbskala, ohne matplotlib
            if pd.isna(v):
                return ""
            t = 0.0 if hi == lo else (v - lo) / (hi - lo)  # 0 = niedriger Cq (dunkel)
            c0, c1 = np.array([13, 54, 107]), np.array([205, 226, 251])
            r, g, b = (c0 + (c1 - c0) * t).astype(int)
            return f"background-color: rgb({r},{g},{b}); color: {'#ffffff' if t < 0.5 else '#0b0b0b'}"

        disp = mat.map(lambda v: "" if pd.isna(v) else f"{v:.2f}")  # Text-Matrix: leere Zellen statt "None"
        st.dataframe(disp.style.apply(lambda _: mat.map(_cell_color), axis=None), use_container_width=True)

# --------------------------------------------------------------------------------------
with tab_label:
    left, right = st.columns([2.3, 1], gap="large")

    with left:
        o1, o2, o3 = st.columns([2.4, 1.1, 1])
        lcolor = o1.radio("Färben nach", ["Gen", "Treatment", "Probe", "Replikat", "Cq"], horizontal=True,
                          key="cb_label")
        st.session_state.sel_mode = o2.radio("Neue Auswahl", ["ersetzen", "hinzufügen"], horizontal=False,
                                             key="sel_mode_radio", label_visibility="collapsed")
        ltext = o3.toggle("Text", value=True, key="txt_label")
        plate_with_selection(d, lcolor, ltext, "plate_label", height=540)
        st.caption("Wells mit **Rechteck-** oder **Lasso-Auswahl** markieren (Werkzeug oben rechts am Plot) "
                   "oder einzeln anklicken. Alternativ rechts einen Bereich eintippen.")

    with right:
        st.subheader("Auswahl")
        rng = st.text_input("Bereich eingeben", placeholder="z.B. A1:H9  ·  A13:P21  ·  B  ·  3:5",
                            help="Mehrere Bereiche mit ; trennen. 'B' = ganze Zeile B, '3:5' = Spalten 3–5, "
                                 "'A:D' = Zeilen A–D.")
        r1, r2 = st.columns(2)
        if r1.button("Übernehmen", use_container_width=True, disabled=not rng):
            try:
                st.session_state.sel = [w for w in parse_range(rng) if w in d.index]
                st.session_state.plot_ver += 1
                st.session_state.last_evt = None
                st.rerun()
            except ValueError as e:
                st.error(str(e))
        if r2.button("Auswahl leeren", use_container_width=True):
            clear_selection()
            st.rerun()
        sel = st.session_state.sel
        st.markdown(f"**Ausgewählt:** {describe_wells(sel)}")

        lab = st.session_state.labels

        def assign(col, value):
            new = st.session_state.labels.copy()
            new.loc[sel, col] = value
            set_labels(new)

        st.markdown("---")
        st.markdown("**① Gen**")
        known_g = sorted({g for g in lab["Gen"] if g})
        g1, g2 = st.columns([2, 1])
        gene = g1.text_input("Gen", key="gene_in", label_visibility="collapsed",
                             placeholder="Gen-Name" + (f" (z.B. {known_g[0]})" if known_g else ""))
        if g2.button("Zuweisen", key="b_gene", use_container_width=True, disabled=not (sel and gene)):
            assign("Gen", gene.strip())
            st.rerun()

        st.markdown("**② Treatment**")
        t1, t2 = st.columns([2, 1])
        treat = t1.text_input("Treatment", key="treat_in", label_visibility="collapsed",
                              placeholder="Treatment-Name")
        if t2.button("Zuweisen", key="b_treat", use_container_width=True, disabled=not (sel and treat)):
            assign("Treatment", treat.strip())
            st.rerun()

        st.markdown("**③ Proben & technische Replikate**")
        p1, p2 = st.columns(2)
        n_rep = p1.number_input("Replikate je Probe", 1, 24, 3)
        start = p2.number_input("Erste Probennummer", 1, 9999, 1)
        rep_dir = st.radio("Replikate liegen …", ["nebeneinander (in einer Zeile)", "untereinander (in einer Spalte)"],
                           horizontal=False)
        order = st.radio("Probennummern laufen zuerst …", ["nach unten ↓", "nach rechts →"], horizontal=True)
        rd = "h" if rep_dir.startswith("neben") else "v"
        od = "down" if order.startswith("nach unten") else "right"
        if sel:
            prev, n_s, inc = sample_blocks(sel, n_rep, rd, od, start)
            first = sorted(prev.items(), key=lambda kv: (kv[1][0], kv[1][1]))[: n_rep * 2]
            ex = " · ".join(f"{w}=P{p}/R{r}" for w, (p, r) in first)
            st.caption(f"→ {n_s} Proben (Nr. {start}–{start + n_s - 1}). Beispiel: {ex} …")
            if inc:
                st.warning(f"{inc} Probe(n) haben weniger als {n_rep} Replikate – Auswahl prüfen.")
        if st.button("Proben nummerieren", type="primary", use_container_width=True, disabled=not sel):
            mapping, _, _ = sample_blocks(sel, n_rep, rd, od, start)
            new = st.session_state.labels.copy()
            for w, (p, r) in mapping.items():
                new.loc[w, "Probe"] = p
                new.loc[w, "Replikat"] = r
            set_labels(new)
            st.rerun()

        st.markdown("---")
        e1, e2 = st.columns(2)
        if e1.button("✕ Ausschließen", use_container_width=True, disabled=not sel):
            assign("Ausschluss", True)
            st.rerun()
        if e2.button("Wieder einschließen", use_container_width=True, disabled=not sel):
            assign("Ausschluss", False)
            st.rerun()
        e3, e4 = st.columns(2)
        if e3.button("🧹 Beschriftung löschen", use_container_width=True, disabled=not sel):
            new = st.session_state.labels.copy()
            new.loc[sel] = empty_labels(sel)
            set_labels(new)
            st.rerun()
        if e4.button("↶ Rückgängig", use_container_width=True,
                     disabled=not st.session_state.get("history")):
            st.session_state.labels = st.session_state.history.pop()
            st.session_state.ver += 1
            st.rerun()

    # Muster kopieren
    with st.expander("📋 Muster auf anderen Plattenbereich kopieren (z.B. gleiches Proben-Layout für ein weiteres Gen)"):
        k1, k2, k3, k4 = st.columns([1.3, 1, 1.6, 1])
        src = k1.text_input("Quellbereich", placeholder="A1:H9", key="cp_src")
        dst = k2.text_input("Ziel: linke obere Ecke", placeholder="A13", key="cp_dst")
        what = k3.multiselect("Was kopieren?", LABEL_COLS[:4], default=["Treatment", "Probe", "Replikat"])
        newgene = k4.text_input("Gen im Ziel (optional)", key="cp_gene")
        if st.button("Kopieren", disabled=not (src and dst and what)):
            try:
                sw = parse_range(src)
                r0 = min(rc(w)[0] for w in sw)
                c0 = min(rc(w)[1] for w in sw)
                dr, dc = _parse_ref(dst)
                if dr is None or dc is None:
                    raise ValueError("Ziel bitte als einzelnes Well angeben, z.B. A13.")
                new = st.session_state.labels.copy()
                n_ok, n_out = 0, 0
                for w in sw:
                    r, c = rc(w)
                    tr, tc = r - r0 + dr, c - c0 + dc
                    if tr >= len(ROWS) or tc > NCOLS:
                        n_out += 1
                        continue
                    tw = f"{ROWS[tr]}{tc}"
                    if tw not in new.index or w not in new.index:
                        continue
                    for col in what:
                        new.loc[tw, col] = st.session_state.labels.loc[w, col]
                    if newgene:
                        new.loc[tw, "Gen"] = newgene.strip()
                    n_ok += 1
                set_labels(new)
                st.success(f"{n_ok} Wells kopiert" + (f", {n_out} lagen außerhalb der Platte." if n_out else "."))
                st.rerun()
            except ValueError as e:
                st.error(str(e))

    # Tabelle zum Nachbearbeiten
    st.subheader("Tabelle (Feinkorrektur)")
    only_sel = st.toggle("Nur ausgewählte Wells zeigen", value=bool(st.session_state.sel))
    tbl = d[["Cq"] + LABEL_COLS].copy()
    if only_sel and st.session_state.sel:
        tbl = tbl.loc[st.session_state.sel]
    edited = st.data_editor(
        tbl, key=f"editor_{st.session_state.ver}_{only_sel}", use_container_width=True, height=360,
        disabled=["Cq"],
        column_config={
            "Cq": st.column_config.NumberColumn(format="%.2f"),
            "Probe": st.column_config.NumberColumn(step=1, min_value=0),
            "Replikat": st.column_config.NumberColumn(step=1, min_value=1),
            "Ausschluss": st.column_config.CheckboxColumn(),
        })
    edited = clean_labels(edited[LABEL_COLS])
    cur = st.session_state.labels.loc[edited.index]
    if not edited.astype(str).equals(cur.astype(str)):
        push_history()
        new = st.session_state.labels.copy()
        new.loc[edited.index] = edited
        st.session_state.labels = clean_labels(new)
        st.rerun()

# --------------------------------------------------------------------------------------
with tab_res:
    summ = replicate_summary(d, sd_thr)
    if summ.empty:
        st.info("Noch keine vollständig beschrifteten Wells (Gen + Treatment + Probe). "
                "Bitte zuerst im Tab **Beschriften** zuordnen.")
        st.stop()

    st.subheader("Technische Replikate")
    n_flag = int((summ["Flag"] != "").sum())
    sug = outlier_suggestions(d, sd_thr)
    f1, f2 = st.columns([3, 1])
    f1.markdown(f"{len(summ)} Gruppen (Gen × Treatment × Probe), davon **{n_flag} auffällig**.")
    if sug and f2.button(f"Ausreißer ausschließen ({len(sug)})", use_container_width=True,
                         help="Pro auffälliger Gruppe (≥3 Replikate) wird das Well mit der größten Abweichung "
                              "vom Median ausgeschlossen – nur wenn die SD danach unter der Schwelle liegt. "
                              f"Betrifft: {', '.join(sug)}"):
        new = st.session_state.labels.copy()
        new.loc[sug, "Ausschluss"] = True
        set_labels(new)
        st.rerun()
    st.dataframe(
        summ.style.format({"Cq_Mittel": "{:.2f}", "Cq_SD": "{:.3f}", "Cq_Min": "{:.2f}", "Cq_Max": "{:.2f}"},
                          na_rep="–")
        .apply(lambda r: ["background-color: rgba(227,73,72,0.15)" if r["Flag"] else "" for _ in r], axis=1),
        use_container_width=True, hide_index=True, height=320)

    st.subheader("Relative Expression (ΔΔCq)")
    per = fcs = None
    genes = sorted(summ["Gen"].unique())
    treats = list(dict.fromkeys(summ["Treatment"]))
    if len(genes) < 2:
        st.info("Für ΔΔCq werden mindestens zwei Gene benötigt (Referenz- und Zielgen).")
    else:
        guess = [g for g in genes if re.search(r"gapdh|actb|b-?actin|18s|rpl|hprt|tbp|ppia|ubc|ywhaz|b2m|36b4",
                                               g, re.I)] or genes[:1]
        q1, q2, q3 = st.columns([2, 1.3, 1])
        refs = q1.multiselect("Referenzgen(e)", genes, default=guess)
        ctrl_guess = next((t for t in treats if re.search(r"ctrl|control|kontroll|veh|dmso|untreat|mock|pbs",
                                                           t, re.I)), treats[0])
        control = q2.selectbox("Kontroll-Treatment", treats, index=treats.index(ctrl_guess))
        log2 = q3.toggle("log2-Skala", value=False)
        if refs and len(refs) < len(genes):
            per, fcs = ddcq(summ, refs, control)
            if per.empty:
                st.warning("Keine passenden Paare aus Ziel- und Referenzgen für dieselbe Probe × Treatment gefunden.")
            else:
                st.plotly_chart(fc_figure(per, fcs, control, log2), use_container_width=True,
                                config={"displaylogo": False})
                st.caption("Balken = Mittelwert der Fold Changes, Fehlerbalken = SD, Punkte = einzelne Proben. "
                           "ΔCq = Cq(Ziel) − Ø Cq(Referenz); ΔΔCq = ΔCq − Ø ΔCq(Kontrolle); FC = 2^−ΔΔCq."
                           + (" p-Wert: Welch-t-Test der ΔCq-Werte gegen die Kontrolle." if _stats else ""))
                st.dataframe(fcs.style.format({c: "{:.3f}" for c in fcs.columns if fcs[c].dtype.kind == "f"},
                                              na_rep="–"), use_container_width=True, hide_index=True)
                with st.expander("Werte je Probe"):
                    st.dataframe(per.style.format({c: "{:.3f}" for c in per.columns if per[c].dtype.kind == "f"},
                                                  na_rep="–"), use_container_width=True, hide_index=True)
        else:
            st.info("Bitte mindestens ein Referenzgen wählen (und mindestens ein Zielgen übrig lassen).")

    st.subheader("Export")
    wells_out = d.reset_index()[["Well", "Cq", "Cq_eff"] + LABEL_COLS].rename(columns={"Cq_eff": "Cq_ausgewertet"})
    st.download_button("⬇️ Ergebnisse als Excel", type="primary",
                       data=to_excel({"Wells": wells_out, "Replikate": summ,
                                      "ddCq_je_Probe": per, "ddCq_Zusammenfassung": fcs}),
                       file_name=f"qPCR_Auswertung_{up.name.rsplit('.', 1)[0]}.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
