from pathlib import Path
import csv

from chembl_webresource_client.new_client import new_client


DATA_DIR = Path(__file__).parent.parent / "data"


def download_chembl_assays_metadata(
    target_chembl_id="CHEMBL274", save_dir=DATA_DIR, filename="ccr5_assays.csv"
):
    activity = new_client.activity
    target_activities = activity.filter(target_chembl_id=target_chembl_id).filter(
        standard_type="IC50"
    )

    columns = target_activities[0].keys()
    csv_file = save_dir / filename
    try:
        with open(csv_file, "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=columns)
            writer.writeheader()
            for data in target_activities:
                writer.writerow(data)
    except IOError as e:
        print(e)
