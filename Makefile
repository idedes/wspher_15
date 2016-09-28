

dist: wspher15.tar.gz

wspher15.tar.gz: Modules15 dens_exper ws_exper ws15_pseudo-exper_IFDENS-0.d \
                  wspher_15.d  wspher_15.f  wspher_15_cmp_ifc  wspher_15_split
	tar zcvf wspher15.tar.gz Modules15 dens_exper ws_exper ws15_pseudo-exper_IFDENS-0.d \
                  wspher_15.d  wspher_15.f  wspher_15_cmp_ifc  wspher_15_split