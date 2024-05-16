from astro_queries.astro_object import AstroFrame
from psql_connector import PgConnector
from unittest.mock import patch


@patch("astro_queries.astro_object.AstroFrame")
def test_astroframe(_):
    model = AstroFrame()
    assert model.detections == None
    assert model.forced_photometry == None
    assert model.detections_filtered == None
    assert model.forced_photometry_filtered == None
    assert isinstance(model, AstroFrame)
    assert isinstance(model.psql_conn, PgConnector)
