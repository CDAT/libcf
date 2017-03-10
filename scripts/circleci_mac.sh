export PATH=${HOME}/miniconda/bin:${PATH}
echo "CIRCLE CI BRANCH:"$CIRCLE_BRANCH
echo "CI_PULL_REQUESTS"$CI_PULL_REQUESTS
echo "CI_PULL_REQUEST"$CI_PULL_REQUEST

export PREFIX=${HOME}/miniconda
export CXXFLAGS="${CXXFLAGS} -fno-common"
export CFLAGS="${XFLAGS} -fno-common"
#export MACOSX_DEPLOYMENT_TARGET=$(sw_vers -productVersion | sed -E "s/([0-9]+\.[0-9]+).*/\1/")
#export DYLD_FALLBACK_LIBRARY_PATH=${PREFIX}/lib
#export DYLD_LIBRARY_PATH=${PREFIX}/lib
export LDFLAGS="${LDFLAGS} -L${PREFIX}/lib -lpython " 
export LFLAGS="-fPIC"
env

pwd
./configure --prefix=${PREFIX}
make
make install

