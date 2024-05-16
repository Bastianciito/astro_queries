import logging


def get_detections_and_force_photometry_by_oid(oid: str = None, **kwargs) -> list[str]:

    if oid is None:
        logging.warning("oid was not provider this must set to get data")
        return None

    astro_query_dd = f"""select * from alerce.detection as d where d.oid = {oid};"""
    astro_query_ff = (
        f"""select * from alerce.forced_photometry as d where d.oid = {oid};"""
    )

    return [astro_query_dd, astro_query_ff]
