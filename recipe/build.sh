export CFLAGS="-Wall -m64 -pipe -O2  -fPIC -I${PREFIX}/include ${CFLAGS}"
export CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
export LFLAGS="-fPIC ${LFLAGS}"
export FC=""

./configure --prefix=${PREFIX}
make
make install
if [ `uname` == Linux ]; then
    LDSHARED="$CC -shared -pthread"  ${PYTHON} setup.py install;
else
    # LDSHARED_FLAGS="-bundle -undefined dynamic_lookup"  ${PYTHON} setup.py install;
    ${PYTHON} setup.py install;
fi
