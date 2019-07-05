#source activate "${CONDA_DEFAULT_ENV}"
export CFLAGS="-Wall -m64 -pipe -O2  -fPIC ${CFLAGS}"
export CXXFLAGS="${CFLAGS} ${CXXFLAGS}"
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS}"
export LDFLAGS="-L${PREFIX}/lib ${LDFLAGS}"
export LFLAGS="-fPIC ${LFLAGS}"
export FC=""
./configure --prefix=${PREFIX}
make
make install
if [ $(uname) == "Linux" ];then
    export LDSHARED="$CC -shared -pthread"
    LDSHARED="$CC -shared -pthread" python setup.py install
else
    if [${HOME} == "/Users/distiller"]; then
        export  CFLAGS="-Wl,-syslibroot / -isysroot / $(CFLAGS)"
    fi
    python setup.py install
fi
