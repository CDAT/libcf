export CFLAGS="-Wall -m64 -pipe -O2  -fPIC -I${PREFIX}/include ${CFLAGS}"
export CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
export LFLAGS="-fPIC ${LFLAGS}"
export FC=""

gcc --version
# needed for clang_osx-64
if [ ${HOME} == "/Users/distiller" ]; then
    echo "xxx xxx xxx gcc --version"
    gcc --version
    export CFLAGS="-Wl,-syslibroot -isysroot ${CFLAGS}"
    # configure need this otherwise "error.h" is not found and configure report netcdf.h 
    export CPPFLAGS="-Wl,-syslibroot -isysroot -I${PREFIX}/include ${CPPFLAGS}"
fi
./configure --prefix=${PREFIX}
make
make install
if [ `uname` == Linux ]; then
    LDSHARED="$CC -shared -pthread"  ${PYTHON} setup.py install;
else
    # LDSHARED_FLAGS="-bundle -undefined dynamic_lookup"  ${PYTHON} setup.py install;
    ${PYTHON} setup.py install;
fi
