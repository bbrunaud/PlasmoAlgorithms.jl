prob = ones(length(scenarios)) / length(scenarios)

	

Sulphur_GO_data = zeros(length(crudes), length(scenarios))
VaccuumResidue_data = zeros(length(crudes), length(scenarios))

Sulphur_GO_data[1,1] = 0.175462361694935
Sulphur_GO_data[2,1] = 0.294652508038040
Sulphur_GO_data[3,1] = 0.162283678416374
Sulphur_GO_data[4,1] = 0.206694887604381
Sulphur_GO_data[5,1] = 0.251898140316648
Sulphur_GO_data[6,1] = 0.580613811870159
Sulphur_GO_data[7,1] = 0.677417307042595
Sulphur_GO_data[8,1] = 1.432271625833148
Sulphur_GO_data[9,1] = 0.330216763244219
Sulphur_GO_data[10,1] = 0.997515844596250

Sulphur_GO_data[1,2] = 0.167549819967135
Sulphur_GO_data[2,2] = 0.290271332261500
Sulphur_GO_data[3,2] = 0.164880392136720
Sulphur_GO_data[4,2] = 0.206586976169398
Sulphur_GO_data[5,2] = 0.253427454865333
Sulphur_GO_data[6,2] = 0.784156500579452
Sulphur_GO_data[7,2] = 0.702454217689116
Sulphur_GO_data[8,2] = 1.662956632895603
Sulphur_GO_data[9,2] = 0.366462109670776
Sulphur_GO_data[10,2] = 1.094567984768590

Sulphur_GO_data[1,3] = 0.164839394038247
Sulphur_GO_data[2,3] = 0.263451555582878
Sulphur_GO_data[3,3] = 0.174766558189615
Sulphur_GO_data[4,3] = 0.180715667985895
Sulphur_GO_data[5,3] = 0.266523539534390
Sulphur_GO_data[6,3] = 0.640118997248147
Sulphur_GO_data[7,3] = 0.821328873482134
Sulphur_GO_data[8,3] = 1.632041974957852
Sulphur_GO_data[9,3] = 0.283336714252063
Sulphur_GO_data[10,3] = 1.110311707777253

Sulphur_GO_data[1,4] = 0.120704871098502
Sulphur_GO_data[2,4] = 0.301789853934891
Sulphur_GO_data[3,4] = 0.151830880633875
Sulphur_GO_data[4,4] = 0.202606738005625
Sulphur_GO_data[5,4] = 0.271531459945764
Sulphur_GO_data[6,4] = 0.634220462682843
Sulphur_GO_data[7,4] = 0.785026256713499
Sulphur_GO_data[8,4] = 1.485454001833434
Sulphur_GO_data[9,4] = 0.361354585155993
Sulphur_GO_data[10,4] = 1.139047001811377

Sulphur_GO_data[1,5] = 0.169209557082088
Sulphur_GO_data[2,5] = 0.320904353761120
Sulphur_GO_data[3,5] = 0.182409949128301
Sulphur_GO_data[4,5] = 0.200657264933307
Sulphur_GO_data[5,5] = 0.261442784852258
Sulphur_GO_data[6,5] = 0.614989673932722
Sulphur_GO_data[7,5] = 0.699533672833922
Sulphur_GO_data[8,5] = 1.584736177999601
Sulphur_GO_data[9,5] = 0.356735777467061
Sulphur_GO_data[10,5] = 1.120244055156392

Sulphur_GO_data[1,6] = 0.167762673197983
Sulphur_GO_data[2,6] = 0.287018324370152
Sulphur_GO_data[3,6] = 0.160773163122308
Sulphur_GO_data[4,6] = 0.199719030590223
Sulphur_GO_data[5,6] = 0.266686826709365
Sulphur_GO_data[6,6] = 0.591783067720569
Sulphur_GO_data[7,6] = 0.721695697299784
Sulphur_GO_data[8,6] = 1.573583896783946
Sulphur_GO_data[9,6] = 0.292898608736301
Sulphur_GO_data[10,6] = 1.220203402690557

Sulphur_GO_data[1,7] = 0.164655479285158
Sulphur_GO_data[2,7] = 0.327476220572340
Sulphur_GO_data[3,7] = 0.159860793611463
Sulphur_GO_data[4,7] = 0.247996625209466
Sulphur_GO_data[5,7] = 0.274299992186458
Sulphur_GO_data[6,7] = 0.736566713864835
Sulphur_GO_data[7,7] = 0.780683050477949
Sulphur_GO_data[8,7] = 1.749874355657192
Sulphur_GO_data[9,7] = 0.301294658686343
Sulphur_GO_data[10,7] = 1.289260370123556

Sulphur_GO_data[1,8] = 0.140104966304307
Sulphur_GO_data[2,8] = 0.216713165757913
Sulphur_GO_data[3,8] = 0.181486631169982
Sulphur_GO_data[4,8] = 0.215582417271255
Sulphur_GO_data[5,8] = 0.233566312984396
Sulphur_GO_data[6,8] = 0.664034994469934
Sulphur_GO_data[7,8] = 0.868796984527227
Sulphur_GO_data[8,8] = 1.607962569310379
Sulphur_GO_data[9,8] = 0.347918915310447
Sulphur_GO_data[10,8] = 1.053006498923610

Sulphur_GO_data[1,9] = 0.148225635593577
Sulphur_GO_data[2,9] = 0.294152054839002
Sulphur_GO_data[3,9] = 0.177212774079495
Sulphur_GO_data[4,9] = 0.176635470167026
Sulphur_GO_data[5,9] = 0.319989228117662
Sulphur_GO_data[6,9] = 0.744359587719000
Sulphur_GO_data[7,9] = 0.788936397907204
Sulphur_GO_data[8,9] = 1.482911297106921
Sulphur_GO_data[9,9] = 0.396627846723051
Sulphur_GO_data[10,9] = 1.105307386048282

Sulphur_GO_data[1,10] = 0.178278616447492
Sulphur_GO_data[2,10] = 0.243548492378170
Sulphur_GO_data[3,10] = 0.182004059334428
Sulphur_GO_data[4,10] = 0.190162521640725
Sulphur_GO_data[5,10] = 0.258513709552509
Sulphur_GO_data[6,10] = 0.640538786184700
Sulphur_GO_data[7,10] = 0.811432512246921
Sulphur_GO_data[8,10] = 1.702688580496033
Sulphur_GO_data[9,10] = 0.367413519283910
Sulphur_GO_data[10,10] = 1.117850008371100


VaccuumResidue_data[1,1] = 0.093106903546301
VaccuumResidue_data[2,1] = 0.253055752661409
VaccuumResidue_data[3,1] = 0.187873253699894
VaccuumResidue_data[4,1] = 0.136104541235664
VaccuumResidue_data[5,1] = 0.263541314032277
VaccuumResidue_data[6,1] = 0.393605035141455
VaccuumResidue_data[7,1] = 0.113254596493256
VaccuumResidue_data[8,1] = 0.250127234251550
VaccuumResidue_data[9,1] = 0.222115836002969
VaccuumResidue_data[10,1] = 0.268252879447202

VaccuumResidue_data[1,2] = 0.097604710038233
VaccuumResidue_data[2,2] = 0.239505438609363
VaccuumResidue_data[3,2] = 0.197364539589748
VaccuumResidue_data[4,2] = 0.117934904612481
VaccuumResidue_data[5,2] = 0.297084692504631
VaccuumResidue_data[6,2] = 0.368591518572553
VaccuumResidue_data[7,2] = 0.102925352969410
VaccuumResidue_data[8,2] = 0.259115310500231
VaccuumResidue_data[9,2] = 0.194907585352156
VaccuumResidue_data[10,2] = 0.267096198950366

VaccuumResidue_data[1,3] = 0.095138912211650
VaccuumResidue_data[2,3] = 0.287071037580188
VaccuumResidue_data[3,3] = 0.180787543551247
VaccuumResidue_data[4,3] = 0.130782922934147
VaccuumResidue_data[5,3] = 0.254282157954429
VaccuumResidue_data[6,3] = 0.401895378801254
VaccuumResidue_data[7,3] = 0.126425673532794
VaccuumResidue_data[8,3] = 0.299619569397980
VaccuumResidue_data[9,3] = 0.195587687830581
VaccuumResidue_data[10,3] = 0.292522317741529

VaccuumResidue_data[1,4] = 0.097061241688740
VaccuumResidue_data[2,4] = 0.284167954353054
VaccuumResidue_data[3,4] = 0.214574188155203
VaccuumResidue_data[4,4] = 0.135041928093711
VaccuumResidue_data[5,4] = 0.220511260300269
VaccuumResidue_data[6,4] = 0.379648534592088
VaccuumResidue_data[7,4] = 0.098442204330274
VaccuumResidue_data[8,4] = 0.279284216838314
VaccuumResidue_data[9,4] = 0.199118639766024
VaccuumResidue_data[10,4] = 0.249136288755436

VaccuumResidue_data[1,5] = 0.107831115165110
VaccuumResidue_data[2,5] = 0.271281063140250
VaccuumResidue_data[3,5] = 0.252246102744745
VaccuumResidue_data[4,5] = 0.148231486771876
VaccuumResidue_data[5,5] = 0.245009012062501
VaccuumResidue_data[6,5] = 0.348552577991371
VaccuumResidue_data[7,5] = 0.113415570264203
VaccuumResidue_data[8,5] = 0.272443644972854
VaccuumResidue_data[9,5] = 0.232821015493706
VaccuumResidue_data[10,5] = 0.275840080171675

VaccuumResidue_data[1,6] = 0.095930323504223
VaccuumResidue_data[2,6] = 0.235613256177243
VaccuumResidue_data[3,6] = 0.199680801550877
VaccuumResidue_data[4,6] = 0.144853410240681
VaccuumResidue_data[5,6] = 0.260372589557835
VaccuumResidue_data[6,6] = 0.364732188023459
VaccuumResidue_data[7,6] = 0.108333348525128
VaccuumResidue_data[8,6] = 0.327933017465411
VaccuumResidue_data[9,6] = 0.201600235598544
VaccuumResidue_data[10,6] = 0.234413312459571

VaccuumResidue_data[1,7] = 0.100973827528024
VaccuumResidue_data[2,7] = 0.227271993296697
VaccuumResidue_data[3,7] = 0.223140327757325
VaccuumResidue_data[4,7] = 0.133492012113300
VaccuumResidue_data[5,7] = 0.247159369536966
VaccuumResidue_data[6,7] = 0.341607063576925
VaccuumResidue_data[7,7] = 0.110153499222355
VaccuumResidue_data[8,7] = 0.220791280338958
VaccuumResidue_data[9,7] = 0.210134935987897
VaccuumResidue_data[10,7] = 0.263135348796899

VaccuumResidue_data[1,8] = 0.095039622146332
VaccuumResidue_data[2,8] = 0.254446089409426
VaccuumResidue_data[3,8] = 0.250346070878934
VaccuumResidue_data[4,8] = 0.102372832191661
VaccuumResidue_data[5,8] = 0.223490374781321
VaccuumResidue_data[6,8] = 0.309993684741986
VaccuumResidue_data[7,8] = 0.113453003895140
VaccuumResidue_data[8,8] = 0.283566386494435
VaccuumResidue_data[9,8] = 0.233476449756469
VaccuumResidue_data[10,8] = 0.224016865054729

VaccuumResidue_data[1,9] = 0.111663955803506
VaccuumResidue_data[2,9] = 0.256882796478658
VaccuumResidue_data[3,9] = 0.205606382356264
VaccuumResidue_data[4,9] = 0.122446424607985
VaccuumResidue_data[5,9] = 0.249067631365377
VaccuumResidue_data[6,9] = 0.337262296403452
VaccuumResidue_data[7,9] = 0.113220077128431
VaccuumResidue_data[8,9] = 0.254136280512678
VaccuumResidue_data[9,9] = 0.175364406163968
VaccuumResidue_data[10,9] = 0.300262218425598

VaccuumResidue_data[1,10] = 0.100788325280492
VaccuumResidue_data[2,10] = 0.302401528667221
VaccuumResidue_data[3,10] = 0.225015907636719
VaccuumResidue_data[4,10] = 0.134127876841188
VaccuumResidue_data[5,10] = 0.243869222905820
VaccuumResidue_data[6,10] = 0.309778853003346
VaccuumResidue_data[7,10] = 0.118286921872142
VaccuumResidue_data[8,10] = 0.260432479499209
VaccuumResidue_data[9,10] = 0.199697212614223
VaccuumResidue_data[10,10] = 0.293026406367485

