from pathlib import Path
import csv

from chembl_webresource_client.new_client import new_client

DATA_DIR = Path("../data")


def get_ccr4_metadata(dir=DATA_DIR):
    activities = new_client.activity
    activities.filter(
        target_chembl_id="CHEMBL2414", pchembl_value__isnull=False, standard_type="IC50"
    )

    columns = activities[0].keys()
    csv_file = DATA_DIR / "ccr4_ic50_meta.csv"
    try:
        with open(csv_file, "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=columns)
            writer.writeheader()
            for data in activities:
                writer.writerow(data)
    except IOError:
        print("I/O error")
