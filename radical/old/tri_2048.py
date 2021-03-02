from src.fp import fp_mul, fp_sqr, fp_exp
def tri_2048(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_101      = fp_mul(_10,_11)
	_111      = fp_mul(_10,_101)
	_1010     = fp_mul(_11,_111)
	_1110     = fp_sqr(_111)
	_1111     = fp_mul(_1,_1110)
	_10000    = fp_mul(_1,_1111)
	_10100    = fp_mul(_101,_1111)
	_11000    = fp_mul(_1010,_1110)
	_11010    = fp_mul(_10,_11000)
	_110010   = fp_mul(_11000,_11010)
	_111001   = fp_mul(_111,_110010)
	_1000011  = fp_mul(_1010,_111001)
	_1000101  = fp_mul(_10,_1000011)
	_1001000  = fp_mul(_11,_1000101)
	_1010001  = fp_mul(_1110,_1000011)
	_1101011  = fp_mul(_11010,_1010001)
	_1111001  = fp_mul(_1110,_1101011)
	_10000000 = fp_mul(_111,_1111001)
	_10000101 = fp_mul(_101,_10000000)
	_11001000 = fp_mul(_1000011,_10000101)
	_11100010 = fp_mul(_11010,_11001000)
	_11101011 = fp_mul(_1101011,_10000000)
	_11101110 = fp_mul(_11,_11101011)
	i26       = fp_mul(_1000011,_11101011)
	i27       = fp_mul(_1111001,i26)
	i28       = fp_mul(_1101011,i27)
	i29       = fp_mul(_111001,i28)
	i30       = fp_mul(_110010,i29)
	i31       = fp_mul(_10000101,i30)
	i32       = fp_mul(_10100,i31)
	i33       = fp_mul(_11101011,i32)
	i34       = fp_mul(_10000000,i33)
	i35       = fp_mul(i29,i34)
	i36       = fp_mul(_11101110,i35)
	x11       = fp_mul(_1000101,i36)
	i38       = fp_mul(_10000,x11)
	i39       = fp_mul(_1010001,i38)
	i40       = fp_mul(_11100010,i39)
	i41       = fp_mul(_1111,i40)
	i42       = fp_mul(_1010,i41)
	i43       = fp_mul(_11001000,i42)
	i44       = fp_mul(i27,i43)
	i45       = fp_mul(i35,i44)
	i46       = fp_mul(i33,i45)
	i47       = fp_mul(_1001000,i46)
	i48       = fp_mul(i26,i47)
	i49       = fp_mul(i42,i48)
	i50       = fp_mul(i43,i49)
	i51       = fp_mul(i31,i50)
	i52       = fp_mul(i38,i51)
	i53       = fp_mul(i44,i52)
	i55       = fp_mul(fp_sqr(i53),i39)
	i56       = fp_mul(i34,i55)
	i57       = fp_mul(i40,i56)
	i58       = fp_mul(i28,i57)
	i59       = fp_mul(i41,i58)
	i60       = fp_mul(i32,i59)
	i61       = fp_mul(i50,i60)
	i62       = fp_mul(i52,i61)
	i63       = fp_mul(i51,i62)
	i64       = fp_mul(i30,i63)
	i65       = fp_mul(i45,i64)
	i66       = fp_mul(i48,i65)
	i67       = fp_mul(i46,i66)
	i68       = fp_mul(i36,i67)
	i69       = fp_mul(i57,i68)
	i70       = fp_mul(i59,i69)
	i71       = fp_mul(i53,i70)
	i72       = fp_mul(i47,i71)
	i73       = fp_mul(i65,i72)
	i74       = fp_mul(i61,i73)
	i75       = fp_mul(i49,i74)
	i76       = fp_mul(i73,i75)
	i77       = fp_mul(i56,i76)
	i78       = fp_mul(i62,i77)
	i79       = fp_mul(i72,i78)
	i80       = fp_mul(i67,i79)
	i81       = fp_mul(i68,i80)
	i82       = fp_mul(i60,i81)
	i83       = fp_mul(i70,i82)
	i84       = fp_mul(i63,i83)
	i85       = fp_mul(i69,i84)
	i86       = fp_mul(i66,i85)
	i87       = fp_mul(i55,i86)
	i88       = fp_mul(i76,i87)
	i89       = fp_mul(i81,i88)
	i90       = fp_mul(i74,i89)
	i91       = fp_mul(i71,i90)
	i92       = fp_mul(i82,i91)
	i93       = fp_mul(i85,i92)
	i94       = fp_mul(i84,i93)
	i95       = fp_mul(i92,i94)
	i96       = fp_mul(i75,i95)
	i97       = fp_mul(i64,i96)
	i98       = fp_mul(i78,i97)
	i99       = fp_mul(i58,i98)
	i100      = fp_mul(i86,i99)
	i101      = fp_mul(i77,i100)
	i102      = fp_mul(i98,i101)
	i103      = fp_mul(i80,i102)
	i104      = fp_mul(i101,i103)
	i105      = fp_mul(i91,i104)
	i106      = fp_mul(i87,i105)
	i107      = fp_mul(i83,i106)
	i108      = fp_mul(i96,i107)
	i109      = fp_mul(i88,i108)
	i110      = fp_mul(i100,i109)
	i111      = fp_mul(i103,i110)
	i112      = fp_mul(i90,i111)
	i113      = fp_mul(i95,i112)
	i114      = fp_mul(i104,i113)
	i115      = fp_mul(i97,i114)
	i116      = fp_mul(i89,i115)
	i117      = fp_mul(i79,i116)
	i118      = fp_mul(i112,i117)
	i119      = fp_mul(i105,i118)
	i120      = fp_mul(i93,i119)
	i121      = fp_mul(i102,i120)
	i123      = fp_mul(fp_sqr(i121),i108)
	i124      = fp_mul(i107,i123)
	i125      = fp_mul(i119,i124)
	i126      = fp_mul(i116,i125)
	i127      = fp_mul(i99,i126)
	i128      = fp_mul(i125,i127)
	i129      = fp_mul(i109,i128)
	i130      = fp_mul(i113,i129)
	i131      = fp_mul(i115,i130)
	i132      = fp_mul(i111,i131)
	i133      = fp_mul(i121,i132)
	i134      = fp_mul(i114,i133)
	i135      = fp_mul(i94,i134)
	i136      = fp_mul(i110,i135)
	x32       = fp_mul(i106,i136)
	i236      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i129,268435456),i117),1),137438953472),i135),1),4294967296)
	i308      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i128,i236),1),34359738368),i120),1),17179869184),i124)
	i410      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i308,34359738368),i133),1),34359738368),i136),1),1073741824)
	i480      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i123,i410),1),17179869184),i126),1),8589934592),i130)
	i583      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i480,17179869184),i132),1),68719476736),i134),1),2147483648)
	i651      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i127,i583),1),1073741824),i118),1),34359738368),i131)
	i749      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i651,4294967296),x32),1),4294967296),x32),1),4294967296)
	i816      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i749),1),4294967296),x32),1),4294967296),x32)
	i914      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i816,4294967296),x32),1),4294967296),x32),1),4294967296)
	i981      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i914),1),4294967296),x32),1),4294967296),x32)
	i1079     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i981,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1146     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1079),1),4294967296),x32),1),4294967296),x32)
	i1244     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1146,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1311     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1244),1),4294967296),x32),1),4294967296),x32)
	i1409     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1311,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1476     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1409),1),4294967296),x32),1),4294967296),x32)
	i1574     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1476,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1641     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1574),1),4294967296),x32),1),4294967296),x32)
	i1739     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1641,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1806     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1739),1),4294967296),x32),1),4294967296),x32)
	i1904     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1806,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1971     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1904),1),4294967296),x32),1),4294967296),x32)
	i2069     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1971,4294967296),x32),1),4294967296),x32),1),4294967296)
	i2136     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i2069),1),4294967296),x32),1),4294967296),x32)
	i2213     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2136,4294967296),x32),1),4294967296),x32),1),2048)
	return       fp_mul(x11,i2213)
