# aquaplanet_ozone_data
Documentation of the aquaplanet ozone data available for CMIP6/CFMIP.

# Low-top version
The original data for the Aqua-Planet Experiment project is available from <http://www.met.reading.ac.uk/~mike/APE/ape_ozone.html>

My low-top version should be equivalent, and is available at <http://dx.doi.org/10.5065/D61834Q6>

## Description
An idealized, time-invariant, zonally averaged, zonally symmetric​ ozone dataset designed by the Aqua-Planet Experiment​ ​​P​​r​​o​​j​​e​​c​​t. Reformatted for use in CAM aqua-planet experiments: The original APE dataset was expanded to the 3 spatial dimensions and "time", mapped to the T42 (128x64) horizontal grid, mapped to 26 levels between 1003.69 and 0.28 hPa, and 12 copies of the time-invariant field created to represent​ a "monthly climatology" within the dataset. ​​​CAM uses the reformatted dataset as a monthly climatological​ ozone boundary dataset. Input is interpolated to the model grid​ and model time as CAM is running.



# High-top extension
Permanent link: <http://dx.doi.org/10.5065/D64X5653>

## Description
-   An idealized, time-invariant, zonally averaged, zonally symmetric ozone dataset based on that used for the Aqua-Planet Experiment Project (<http://www.met.reading.ac.uk/~mike/APE/ape_spec.html>, see also doi:10.5065/D61834Q6), but expanded vertically for use in high-topped models.

-   To allow high-top models to utilize the dataset, a WACCM ozone data set was blended with the original APE ozone data set. The WACCM data (`/glade/p/cesm/chwg_dev/tilmes/climo/oxid_ozone_WACCM_CCMI_1850_2010_serial_3D_monthly.nc`, provided by Simone Tilmes) is from an Earth configuration. A long-term climatology is constructed by averaging in time and across longitudes, and made symmetric about the equator. The original APE data set and this WACCM climatology are blended by interpolating the APE data set to the WACCM grid and using the interpolated APE data to mask the WACCM data. This effectively uses the WACCM data for upper levels and retains the APE data at the lower levels. The data is then converted from hybrid-sigma levels to pressure levels and expanded to cover 3 spatial dimensions and "time."

-   The data set is mapped to the CAM finite volume nominal 2-degree (144x96) horizontal grid
-   The vertical grid is 69 pressure levels between 1000.00 and 5e-06 hPa
-   12 copies of the time-invariant field created to represent a "monthly climatology" within the dataset
-   CAM uses the reformatted dataset as a monthly climatological ozone boundary dataset. Input is interpolated to the model grid and model time as CAM is running.

## Scripts
The script `new_aqua_o3.ncl` produced the original version, but neglected to save the zonal average. This *has not* been corrected in this repository.

The jupyter notebook `Aquaplanet_Ozone_Revision.ipynb` takes the old version, applies the zonal average, updates the data set, and saves it. It shows a couple of plots just to show that the operation worked.

## History
Version 2 was published on 4/12/2018. The previous version (`aquaplanet_ozone_hightop_c160920.nc​`) was found to be in error and has been replaced with `aquaplanet_ozone_hightop_c20180412.nc`. Thanks to Yoko Tsushima for bringing the error to my attention. Thanks to Gary Strand for assistance in publishing the revised version to Earth System Grid.
