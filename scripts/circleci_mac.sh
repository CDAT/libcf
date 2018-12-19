export PATH=${HOME}/miniconda/bin:${PATH}
echo "CIRCLE CI BRANCH:"$CIRCLE_BRANCH
echo "CI_PULL_REQUESTS"$CI_PULL_REQUESTS
echo "CI_PULL_REQUEST"$CI_PULL_REQUEST

export PREFIX=${HOME}/miniconda
export CXXFLAGS="${CXXFLAGS} -fno-common"
export CFLAGS="${XFLAGS} -fno-common"
# somehow looks for fortran but doesn't seem to need it it's pure C
# looking for fortran makes mac fail
export FC="" 
env

pwd
./configure --prefix=${PREFIX}
make
make install

