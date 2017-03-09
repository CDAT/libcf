#
# This script copies the grid to four files which will have bogus data added
# to them.
#
# $Id: $
#

if [ -f "tst_two_tiles_add_id_bog00.nc" ] 
then
  rm tst_two_tiles_add_id_bog00.nc
fi
if [ -f "tst_two_tiles_add_id_bog01.nc" ] 
then
rm tst_two_tiles_add_id_bog01.nc
fi
if [ -f "tst_two_tiles_add_id_bog10.nc" ] 
then
  rm tst_two_tiles_add_id_bog10.nc
fi
if [ -f "tst_two_tiles_add_id_bog11.nc" ] 
then
  rm tst_two_tiles_add_id_bog11.nc
fi

cp tst_two_tiles_add_id_grid0.nc tst_two_tiles_add_id_bog00.nc
cp tst_two_tiles_add_id_grid0.nc tst_two_tiles_add_id_bog01.nc
cp tst_two_tiles_add_id_grid1.nc tst_two_tiles_add_id_bog10.nc
cp tst_two_tiles_add_id_grid1.nc tst_two_tiles_add_id_bog11.nc
