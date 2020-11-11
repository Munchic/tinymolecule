from pathlib import Path
import csv

from chembl_webresource_client.new_client import new_client

DATA_DIR = Path("/Users/Munchic/Developer/Capstone/tinymolecule/data")


def get_target_metadata(target_chembl_id="CHEMBL274", dir=DATA_DIR):
    activity = new_client.activity
    target_activities = activity.filter(target_chembl_id=target_chembl_id).filter(
        standard_type="IC50"
    )

    columns = target_activities[0].keys()
    csv_file = DATA_DIR / "ccr5_ic50_meta.csv"
    try:
        with open(csv_file, "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=columns)
            writer.writeheader()
            for data in target_activities:
                writer.writerow(data)
    except IOError as e:
        print(e)


get_target_metadata()