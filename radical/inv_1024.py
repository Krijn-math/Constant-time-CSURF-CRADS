from src.fp import fp_mul, fp_sqr, fp_exp
def inv_1024(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_101      = fp_mul(_10,_11)
	_1000     = fp_mul(_11,_101)
	_1011     = fp_mul(_11,_1000)
	_1101     = fp_mul(_10,_1011)
	_1111     = fp_mul(_10,_1101)
	_10000    = fp_mul(_1,_1111)
	_10100    = fp_mul(_101,_1111)
	_10110    = fp_mul(_10,_10100)
	_11000    = fp_mul(_10,_10110)
	_11111    = fp_mul(_1011,_10100)
	_101110   = fp_mul(_1111,_11111)
	_111011   = fp_mul(_1101,_101110)
	_1000000  = fp_mul(_101,_111011)
	_1001000  = fp_mul(_1000,_1000000)
	_1001010  = fp_mul(_10,_1001000)
	_1010000  = fp_mul(_1000,_1001000)
	_1010011  = fp_mul(_11,_1010000)
	_1011000  = fp_mul(_101,_1010011)
	_10000001 = fp_mul(_101110,_1010011)
	_10100000 = fp_mul(_11111,_10000001)
	_10101000 = fp_mul(_1000,_10100000)
	_10110110 = fp_mul(_10110,_10100000)
	_10110111 = fp_mul(_1,_10110110)
	_10111010 = fp_mul(_11,_10110111)
	_11000000 = fp_mul(_11000,_10101000)
	_11101110 = fp_mul(_101110,_11000000)
	i29       = fp_mul(_10110,_11101110)
	i30       = fp_mul(_10000,i29)
	i31       = fp_mul(_1000000,_11101110)
	i32       = fp_sqr(_10100000)
	i33       = fp_mul(_10101000,_10111010)
	i34       = fp_mul(_10111010,i33)
	i35       = fp_mul(_10100,i34)
	i36       = fp_mul(_1001010,i35)
	i37       = fp_mul(_1001010,i36)
	i38       = fp_mul(_1011,i37)
	i39       = fp_mul(_10000001,i38)
	i40       = fp_mul(i32,i34)
	i41       = fp_mul(_1011000,i40)
	i42       = fp_mul(_11000000,i41)
	i43       = fp_mul(_1001000,i42)
	i44       = fp_mul(i29,i43)
	i45       = fp_mul(_10110110,i44)
	i46       = fp_mul(_10,i45)
	i47       = fp_mul(_10000,i45)
	i48       = fp_mul(_11000,i47)
	i49       = fp_mul(_10101000,i48)
	i50       = fp_mul(_10110111,i49)
	i51       = fp_mul(_111011,i50)
	i52       = fp_mul(_1111,i51)
	i53       = fp_mul(_1010011,i52)
	i54       = fp_mul(_11101110,i53)
	i55       = fp_mul(i40,i48)
	i56       = fp_mul(i33,i55)
	i57       = fp_mul(i39,i56)
	i58       = fp_mul(_11,i57)
	i59       = fp_mul(_1001000,i57)
	i60       = fp_mul(i51,i59)
	i61       = fp_mul(i34,i60)
	i62       = fp_mul(i38,i61)
	i63       = fp_mul(_1111,i62)
	i64       = fp_mul(i58,i63)
	i65       = fp_mul(i31,i64)
	i66       = fp_mul(i56,i65)
	i67       = fp_mul(i59,i66)
	i68       = fp_mul(i40,i67)
	i69       = fp_mul(i55,i68)
	i70       = fp_mul(i44,i69)
	i71       = fp_mul(i41,i70)
	i72       = fp_mul(_1000000,i71)
	i73       = fp_mul(i63,i72)
	i74       = fp_mul(i47,i73)
	i75       = fp_mul(i54,i74)
	i76       = fp_mul(i42,i75)
	i77       = fp_mul(i35,i76)
	i78       = fp_mul(i45,i77)
	i79       = fp_mul(i57,i78)
	i80       = fp_sqr(i69)
	i81       = fp_mul(i43,i79)
	i82       = fp_mul(i32,i81)
	i83       = fp_mul(i53,i82)
	i84       = fp_mul(i61,i83)
	i85       = fp_mul(i60,i84)
	i86       = fp_mul(i37,i85)
	i87       = fp_mul(i36,i86)
	i88       = fp_mul(i46,i87)
	i89       = fp_mul(i30,i88)
	i90       = fp_mul(_10100000,i89)
	i91       = fp_mul(i49,i90)
	i92       = fp_mul(i48,i91)
	x16       = fp_mul(_1010000,i92)
	i149      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i80,262144),i76),1),1048576),i78),1),65536)
	i186      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i90,i149),1),65536),i84),1),262144),i89)
	i233      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i186,131072),i75),1),65536),i91),1),4096)
	i273      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i52,i233),1),262144),i66),1),524288),i77)
	i327      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i273,262144),i82),1),65536),i74),1),262144)
	i361      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i79,i327),1),65536),i83),1),32768),i68)
	i421      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i361,4194304),i92),1),524288),i81),1),131072)
	i457      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i88,i421),1),65536),i64),1),131072),i71)
	i511      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i457,1048576),i67),1),131072),i86),1),32768)
	i548      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i72,i511),1),131072),i73),1),131072),i85)
	i600      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i548,32768),i65),1),262144),i70),1),131072)
	i639      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i62,i600),1),1048576),i87),1),65536),x16)
	i689      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i639,65536),x16),1),65536),x16),1),65536)
	i724      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x16,i689),1),65536),x16),1),65536),x16)
	i774      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i724,65536),x16),1),65536),x16),1),65536)
	i809      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x16,i774),1),65536),x16),1),65536),x16)
	i859      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i809,65536),x16),1),65536),x16),1),65536)
	i894      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x16,i859),1),65536),x16),1),65536),x16)
	i944      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i894,65536),x16),1),65536),x16),1),65536)
	i979      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x16,i944),1),65536),x16),1),65536),x16)
	i1029     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i979,65536),x16),1),65536),x16),1),65536)
	i1064     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x16,i1029),1),65536),x16),1),65536),x16)
	i1114     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1064,65536),x16),1),65536),x16),1),65536)
	i1149     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x16,i1114),1),65536),x16),1),65536),x16)
	return       fp_mul(fp_exp(i1149,2048),i50)
