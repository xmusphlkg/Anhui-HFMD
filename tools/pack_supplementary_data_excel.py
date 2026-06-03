from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from pypinyin import lazy_pinyin


ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = ROOT / "outcome" / "supplementary_tables"
OUTPUT_FILE = ROOT / "outcome" / "supplementary_data.xlsx"


DATA_GROUPS = [
    {
        "data_id": "D1",
        "title": "yearly_typed_serotype_composition",
        "sheets": [
            ("typed_composition_yearly", "typed_composition_yearly.csv"),
        ],
    },
    {
        "data_id": "D2",
        "title": "extended_diagnostics_nb_its",
        "sheets": [
            ("its_nb_diag_ext", "its_nb_diagnostics_extended.csv"),
        ],
    },
    {
        "data_id": "D3",
        "title": "typing_fraction_sensitivity",
        "sheets": [
            ("its_typ_frac_sens", "its_typing_fraction_sensitivity.csv"),
        ],
    },
    {
        "data_id": "D4",
        "title": "glm_apc_sensitivity",
        "sheets": [
            ("glm_apc_sens", "joinpoint_glm_apc_sensitivity.csv"),
        ],
    },
    {
        "data_id": "D5",
        "title": "severe_mild_ev71",
        "sheets": [
            ("severe_mild_ev71", "severe_mild_ev71_yearly.csv"),
            ("severe_mild_or", "severe_vs_mild_ev71_or_yearly.csv"),
        ],
    },
    {
        "data_id": "D6",
        "title": "sampling_bias_logit",
        "sheets": [
            ("sampling_bias_logit", "severe_sampling_bias_logit.csv"),
        ],
    },
    {
        "data_id": "D7",
        "title": "county_typed_confidence",
        "sheets": [
            ("county_typed_conf", "supplementary_data_D7_county_map_masking.csv"),
        ],
    },
    {
        "data_id": "D8",
        "title": "satscan_sensitivity",
        "sheets": [
            ("satscan_sens_sum", "satscan_cluster_size_sensitivity_summary.csv"),
            ("satscan_sens_det", "satscan_cluster_size_sensitivity_detail.csv"),
        ],
    },
    {
        "data_id": "D9",
        "title": "postpandemic_age_shift",
        "sheets": [
            ("annual_age_comp", "supplementary_data_D9_annual_age_composition.csv"),
            ("age_shift_period", "supplementary_data_D9_period_summary.csv"),
            ("age_shift_stats", "supplementary_data_D9_statistics.csv"),
        ],
    },
    {
        "data_id": "D10",
        "title": "severe_death_cfr",
        "sheets": [
            ("severe_death_yearly", "severe_death_cfr_yearly.csv"),
            ("severe_death_overall", "severe_death_cfr_overall.csv"),
        ],
    },
]


COUNTY_NAME_OVERRIDES = {
    "和县": "He County",
    "寿县": "Shou County",
    "泗县": "Si County",
    "萧县": "Xiao County",
    "歙县": "She County",
    "黟县": "Yi County",
    "枞阳县": "Zongyang County",
    "埇桥区": "Yongqiao District",
    "谯城区": "Qiaocheng District",
    "涡阳县": "Guoyang County",
    "颍上县": "Yingshang County",
    "颍东区": "Yingdong District",
    "颍州区": "Yingzhou District",
    "颍泉区": "Yingquan District",
    "郊区": "Jiao District",
}


def truncate_sheet_name(name: str) -> str:
    return name[:31]


def to_english_place_name(name: str) -> str:
    if name in COUNTY_NAME_OVERRIDES:
        return COUNTY_NAME_OVERRIDES[name]

    suffix_map = {
        "市": "City",
        "县": "County",
        "区": "District",
    }

    suffix = name[-1] if name and name[-1] in suffix_map else ""
    stem = name[:-1] if suffix else name
    stem_pinyin = "".join(lazy_pinyin(stem)).capitalize()

    if suffix:
        return f"{stem_pinyin} {suffix_map[suffix]}"
    return stem_pinyin


def transform_frame(group: dict[str, object], frame: pd.DataFrame) -> pd.DataFrame:
    if str(group["data_id"]) == "D7" and "name" in frame.columns:
        frame = frame.copy()
        frame["name"] = frame["name"].astype(str).map(to_english_place_name)
    return frame


def autofit_columns(workbook_path: Path) -> None:
    workbook = load_workbook(workbook_path)
    try:
        for worksheet in workbook.worksheets:
            for column_cells in worksheet.columns:
                values = ["" if cell.value is None else str(cell.value) for cell in column_cells]
                max_length = max((len(value) for value in values), default=0)
                column_letter = column_cells[0].column_letter
                worksheet.column_dimensions[column_letter].width = min(max_length + 2, 80)
        workbook.save(workbook_path)
    finally:
        workbook.close()


def write_group_sheet(writer: pd.ExcelWriter, group: dict[str, object]) -> None:
    sheet_name = truncate_sheet_name(str(group["data_id"]))
    start_row = 0

    for table_name, csv_name in group["sheets"]:
        csv_path = SOURCE_DIR / csv_name
        if not csv_path.exists():
            raise FileNotFoundError(f"Missing source CSV: {csv_path}")

        frame = pd.read_csv(csv_path)
        frame = transform_frame(group, frame)

        # Put multiple tables from the same supplementary data block
        # into one worksheet, separated by a short title row and spacing.
        pd.DataFrame([[table_name]]).to_excel(
            writer,
            sheet_name=sheet_name,
            startrow=start_row,
            index=False,
            header=False,
        )
        frame.to_excel(
            writer,
            sheet_name=sheet_name,
            startrow=start_row + 1,
            index=False,
        )
        start_row += len(frame) + 4


def main() -> None:
    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
        for group in DATA_GROUPS:
            write_group_sheet(writer, group)

    autofit_columns(OUTPUT_FILE)
    print(f"Created workbook: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
