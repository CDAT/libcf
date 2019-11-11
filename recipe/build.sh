export CFLAGS="-Wall -m64 -pipe -O2  -fPIC -I${PREFIX}/include ${CFLAGS}"
export CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
export CPPFLAGS="-I${PREFIX}/include -I/usr/local/netcdf/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
export LFLAGS="-fPIC ${LFLAGS}"
export FC=""

echo "CPPFLAGS: $CPPFLAGS"
# needed for clang_osx-64
if [ ${HOME} == "/Users/distiller" ]; then
    export CFLAGS="-Wl,-syslibroot / -isysroot / ${CFLAGS}"
    # configure need this otherwise "error.h" is not found and configure report netcdf.h 
    export CPPFLAGS="-Wl,-syslibroot / -isysroot / -I${PREFIX}/include ${CPPFLAGS}"
fi
./configure --prefix=${PREFIX}
make
make install
if [ $(uname) == "Linux" ];then
    export LDSHARED="$CC -shared -pthread"
    LDSHARED="$CC -shared -pthread" python setup.py install
else
    python setup.py install
fi
