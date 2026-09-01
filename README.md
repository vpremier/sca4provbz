# sca4provbz

Tools for producing and analysing snow-cover products for South Tyrol for the Risk Monitoring project.

## Main Repositories

- `viirs/`: downloads NASA VNP10A1F version 2 granules and converts them
  into daily snow-cover-fraction (SCF) GeoTIFFs at 250 m resolution. The
  output grid uses EPSG:32632 and the extent required by the Province of
  Bolzano. Water pixels are assigned the unavailable value `205` using the
  pre-aligned mask `aux/Water_Mask_aligned.tif`.
- `snow_bullettin/`: updates the daily snow-cover-area (SCA) time series,
  creates the operational snow bulletin, and computes monthly, quarterly and
  yearly snow-cover-duration (SCD) products.
- `snow_bullettin_old/` and `viirs_old/`: previous workflow versions retained
  for reference.


