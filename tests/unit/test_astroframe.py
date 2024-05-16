from astro_queries.astro_object import AstroFrame


def test_astroframe():
    model = AstroFrame()
    assert model.detections == None
    assert model.forced_photometry == None
    assert model.detections_filtered == None
    assert model.forced_photometry_filtered == None
