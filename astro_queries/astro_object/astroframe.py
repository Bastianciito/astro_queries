from ..queries import get_detections_and_force_photometry_by_oid
from psql_connector import PgConnector
import matplotlib.pyplot as plt


class AstroFrame:

    def __init__(self, oid: str = None, **kwargs) -> None:

        self.detections = None
        self.forced_photometry = None
        self.detections_filtered = None
        self.forced_photometry_filtered = None

    def deduplicate(self):

        self.detections_filtered = None
        self.forced_photometry_filtered = None

    def plot_curves(self):

        pass
