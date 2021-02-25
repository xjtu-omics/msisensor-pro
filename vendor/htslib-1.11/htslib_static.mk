HTSLIB_static_LDFLAGS = -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/pengjia/miniconda3/lib -Wl,-rpath-link,/home/pengjia/miniconda3/lib -L/home/pengjia/miniconda3/lib
HTSLIB_static_LIBS = -lpthread -lz -lm -lbz2 -llzma -ldeflate -lcurl -lcrypto
