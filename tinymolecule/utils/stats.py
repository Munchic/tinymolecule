import numpy as np
from scipy import stats


def stats_report(baff_data):
    report = {}
    report["mean"] = np.mean(baff_data)
    report["median"] = np.mean(baff_data)
    report["std"] = np.std(baff_data)
    report["var"] = report["std"] ** 2
    report["95ci"] = stats.norm.interval(0.95, loc=report["mean"], scale=report["std"])

    return report