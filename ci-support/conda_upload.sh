#!/usr/bin/env bash
PKG_NAME=libcf
USER=cdat
export VERSION="3.1.0"
echo "Trying to upload to conda"
echo ""
echo "Activating base env"
source activate base
echo "Making sure conda-build is installed"
conda install "conda-build"
echo "Updating conda"
conda update -y -q conda
if [ `uname` == "Linux" ]; then
    OS=linux-64
    echo "Linux OS"
else
    echo "Mac OS"
    OS=osx-64
fi

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=${HOME}/conda-bld
echo "Cloning recipes"
git clone git://github.com/CDAT/conda-recipes
cd conda-recipes
# uvcdat creates issues for build -c uvcdat confises package and channel
if [[ -d uvcdat ]]; then
    rm -rf uvcdat
fi
if [[ -d libcf ]]; then
    rm -rf libcf 
fi
ln -s ../recipe libcf
export BRANCH=${CIRCLE_BRANCH}
python ./prep_for_build.py  -b ${BRANCH}
if [[ $PY_VER = 'py2' ]]; then
    conda build $PKG_NAME -c conda-forge -c cdat --python=2.7
else
    conda build $PKG_NAME -c conda-forge -c cdat --python=3.6
fi
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l $LABEL $CONDA_BLD_PATH/$OS/$PKG_NAME-$VERSION.`date +%Y*`0.tar.bz2 --force
