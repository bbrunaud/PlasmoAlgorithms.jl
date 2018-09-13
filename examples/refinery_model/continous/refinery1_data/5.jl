prob = ones(length(scenarios)) / length(scenarios)

	

Sulphur_GO_data = zeros(length(crudes), length(scenarios))
VaccuumResidue_data = zeros(length(crudes), length(scenarios))

Sulphur_GO_data[1,1] = 0.142088513077039
Sulphur_GO_data[2,1] = 0.310719599235616
Sulphur_GO_data[3,1] = 0.166528078567320
Sulphur_GO_data[4,1] = 0.212905808207417
Sulphur_GO_data[5,1] = 0.266482491043895
Sulphur_GO_data[6,1] = 0.761729783608670
Sulphur_GO_data[7,1] = 0.634362056891616
Sulphur_GO_data[8,1] = 1.576735109424405
Sulphur_GO_data[9,1] = 0.291490315101932
Sulphur_GO_data[10,1] = 1.089348038245218

Sulphur_GO_data[1,2] = 0.141118541282429
Sulphur_GO_data[2,2] = 0.281477552818021
Sulphur_GO_data[3,2] = 0.167414406429630
Sulphur_GO_data[4,2] = 0.184276995676455
Sulphur_GO_data[5,2] = 0.263494796995319
Sulphur_GO_data[6,2] = 0.692304963412679
Sulphur_GO_data[7,2] = 0.773320932383037
Sulphur_GO_data[8,2] = 1.309305562100601
Sulphur_GO_data[9,2] = 0.306919506227801
Sulphur_GO_data[10,2] = 1.151557551730689

Sulphur_GO_data[1,3] = 0.187121509799093
Sulphur_GO_data[2,3] = 0.306887871469947
Sulphur_GO_data[3,3] = 0.147750555914438
Sulphur_GO_data[4,3] = 0.216077618904068
Sulphur_GO_data[5,3] = 0.253660826161628
Sulphur_GO_data[6,3] = 0.685002013445849
Sulphur_GO_data[7,3] = 0.687156971268290
Sulphur_GO_data[8,3] = 1.538931836234552
Sulphur_GO_data[9,3] = 0.324073104940641
Sulphur_GO_data[10,3] = 1.100321154941514

Sulphur_GO_data[1,4] = 0.138170604978298
Sulphur_GO_data[2,4] = 0.324808881357782
Sulphur_GO_data[3,4] = 0.132956277603478
Sulphur_GO_data[4,4] = 0.194502792675285
Sulphur_GO_data[5,4] = 0.228934729989087
Sulphur_GO_data[6,4] = 0.745981922602771
Sulphur_GO_data[7,4] = 0.930200382982935
Sulphur_GO_data[8,4] = 1.766616103121546
Sulphur_GO_data[9,4] = 0.333366587927294
Sulphur_GO_data[10,4] = 1.155344266991180

Sulphur_GO_data[1,5] = 0.160119037177764
Sulphur_GO_data[2,5] = 0.237070843250389
Sulphur_GO_data[3,5] = 0.179241534916566
Sulphur_GO_data[4,5] = 0.198547381282052
Sulphur_GO_data[5,5] = 0.274516556420276
Sulphur_GO_data[6,5] = 0.631873427540527
Sulphur_GO_data[7,5] = 0.806307527204965
Sulphur_GO_data[8,5] = 1.625188062179657
Sulphur_GO_data[9,5] = 0.369783597954369
Sulphur_GO_data[10,5] = 1.144253677136521


VaccuumResidue_data[1,1] = 0.098476746897009
VaccuumResidue_data[2,1] = 0.254841406932241
VaccuumResidue_data[3,1] = 0.173170328760428
VaccuumResidue_data[4,1] = 0.141239119577067
VaccuumResidue_data[5,1] = 0.255427719783213
VaccuumResidue_data[6,1] = 0.441044400649980
VaccuumResidue_data[7,1] = 0.112523799680884
VaccuumResidue_data[8,1] = 0.276969737693296
VaccuumResidue_data[9,1] = 0.189320400668474
VaccuumResidue_data[10,1] = 0.226116053589760

VaccuumResidue_data[1,2] = 0.117914428890899
VaccuumResidue_data[2,2] = 0.275227182331483
VaccuumResidue_data[3,2] = 0.236554126564275
VaccuumResidue_data[4,2] = 0.124496284282309
VaccuumResidue_data[5,2] = 0.253628478519632
VaccuumResidue_data[6,2] = 0.414977511926427
VaccuumResidue_data[7,2] = 0.096313256035774
VaccuumResidue_data[8,2] = 0.289686439323075
VaccuumResidue_data[9,2] = 0.205015116160727
VaccuumResidue_data[10,2] = 0.264009751721092

VaccuumResidue_data[1,3] = 0.094777581750325
VaccuumResidue_data[2,3] = 0.295376717751326
VaccuumResidue_data[3,3] = 0.177664240429431
VaccuumResidue_data[4,3] = 0.131236523229636
VaccuumResidue_data[5,3] = 0.228221555314528
VaccuumResidue_data[6,3] = 0.349226110581922
VaccuumResidue_data[7,3] = 0.103454097062942
VaccuumResidue_data[8,3] = 0.205819359076627
VaccuumResidue_data[9,3] = 0.177501118103415
VaccuumResidue_data[10,3] = 0.280661872199777

VaccuumResidue_data[1,4] = 0.092419309390921
VaccuumResidue_data[2,4] = 0.272437946676044
VaccuumResidue_data[3,4] = 0.221951668712980
VaccuumResidue_data[4,4] = 0.111483940457563
VaccuumResidue_data[5,4] = 0.242592773579144
VaccuumResidue_data[6,4] = 0.341326361986696
VaccuumResidue_data[7,4] = 0.124169710255900
VaccuumResidue_data[8,4] = 0.287409969852358
VaccuumResidue_data[9,4] = 0.206129072970268
VaccuumResidue_data[10,4] = 0.245554637551950

VaccuumResidue_data[1,5] = 0.118213633855153
VaccuumResidue_data[2,5] = 0.299551972362780
VaccuumResidue_data[3,5] = 0.220937083136309
VaccuumResidue_data[4,5] = 0.144708783532668
VaccuumResidue_data[5,5] = 0.237201293280314
VaccuumResidue_data[6,5] = 0.376551360490621
VaccuumResidue_data[7,5] = 0.104329719484758
VaccuumResidue_data[8,5] = 0.256205294490090
VaccuumResidue_data[9,5] = 0.193152661395673
VaccuumResidue_data[10,5] = 0.273551035003146

