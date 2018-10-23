#source activate "${CONDA_DEFAULT_ENV}"
export CFLAGS="-Wall -m64 -pipe -O2  -fPIC ${CFLAGS}"
export CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
export LFLAGS="-fPIC ${LFLAGS}"
export FC=""
export LDSHARED="$CC -shared -pthread" 
./configure --prefix=${PREFIX}
make
make install
if [ `uname` == Linux ]; then
    LDSHARED="$CC -shared -pthread"  ${PYTHON} setup.py install;
else
    LDSHARED_FLAGS="-bundle -undefined dynamic_lookup"  ${PYTHON} setup.py install;
fi
# if [ `uname` == Darwin ]; then install_name_tool -change /System/Library/Frameworks/Python.framework/Versions/2.7/Python @rpath/libpython2.7.dylib ${SP_DIR}/pycf/*.so ; fi
