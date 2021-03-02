from src.fp import fp_mul, fp_sqr, fp_exp
def inv_1792(x):
	_1 = x
	_10       = fp_sqr(_1)
	_100      = fp_sqr(_10)
	_101      = fp_mul(_1,_100)
	_110      = fp_mul(_1,_101)
	_1000     = fp_mul(_10,_110)
	_1010     = fp_mul(_10,_1000)
	_1101     = fp_mul(_101,_1000)
	_10000    = fp_mul(_110,_1010)
	_10111    = fp_mul(_1010,_1101)
	_11010    = fp_mul(_1010,_10000)
	_100111   = fp_mul(_1101,_11010)
	_101001   = fp_mul(_10,_100111)
	_101100   = fp_mul(_101,_100111)
	_110011   = fp_mul(_1010,_101001)
	_111000   = fp_mul(_101,_110011)
	_111100   = fp_mul(_100,_111000)
	_1011111  = fp_mul(_100111,_111000)
	_10000110 = fp_mul(_100111,_1011111)
	_10110010 = fp_mul(_101100,_10000110)
	_10111111 = fp_mul(_1101,_10110010)
	i21       = fp_mul(_1011111,_10111111)
	i22       = fp_mul(_10000,i21)
	i23       = fp_mul(_101001,i22)
	i24       = fp_mul(_110,i23)
	i25       = fp_mul(_1010,i24)
	i26       = fp_mul(_110,i25)
	i27       = fp_mul(_1011111,i23)
	i28       = fp_mul(_111000,i27)
	i29       = fp_mul(_111100,i28)
	i30       = fp_mul(_11010,i29)
	i31       = fp_mul(_10111,i30)
	i32       = fp_mul(i22,i31)
	i33       = fp_mul(_110,i32)
	i34       = fp_mul(_111000,i33)
	i35       = fp_mul(i25,i34)
	i36       = fp_mul(i21,i35)
	i37       = fp_mul(i27,i36)
	i38       = fp_mul(_110011,i37)
	i39       = fp_mul(i31,i38)
	i40       = fp_mul(_10000110,i39)
	i41       = fp_mul(_1000,i40)
	i42       = fp_mul(_10110010,i41)
	i43       = fp_mul(i29,i42)
	i44       = fp_mul(_10111111,i43)
	i45       = fp_mul(i34,i44)
	i46       = fp_mul(i28,i45)
	i47       = fp_mul(i23,i46)
	i48       = fp_mul(i36,i47)
	i49       = fp_mul(i32,i48)
	i50       = fp_mul(i37,i49)
	i51       = fp_mul(i47,i50)
	i52       = fp_mul(i38,i51)
	i53       = fp_mul(i39,i52)
	i54       = fp_mul(i44,i53)
	i55       = fp_mul(i45,i54)
	i56       = fp_mul(i33,i55)
	i57       = fp_mul(i41,i56)
	i58       = fp_mul(i35,i57)
	i59       = fp_mul(i46,i58)
	i60       = fp_mul(i30,i59)
	i61       = fp_mul(i56,i60)
	i62       = fp_mul(i42,i61)
	i63       = fp_mul(i53,i62)
	i64       = fp_mul(i26,i63)
	i65       = fp_mul(i52,i64)
	i66       = fp_mul(i54,i65)
	i67       = fp_mul(i57,i66)
	i68       = fp_mul(i49,i67)
	i69       = fp_mul(i51,i68)
	i70       = fp_mul(i24,i69)
	i71       = fp_mul(i50,i70)
	i72       = fp_mul(i62,i71)
	i73       = fp_mul(i48,i72)
	i74       = fp_mul(i55,i73)
	i75       = fp_mul(i40,i74)
	i76       = fp_mul(i59,i75)
	i77       = fp_mul(i60,i76)
	i78       = fp_mul(i61,i77)
	i79       = fp_mul(i63,i78)
	i80       = fp_mul(i71,i79)
	i81       = fp_mul(i68,i80)
	i82       = fp_mul(i65,i81)
	i83       = fp_mul(i58,i82)
	i84       = fp_mul(i78,i83)
	i85       = fp_mul(i69,i84)
	i86       = fp_mul(i72,i85)
	i87       = fp_mul(i67,i86)
	i89       = fp_mul(fp_sqr(i87),i66)
	i90       = fp_mul(i80,i89)
	i91       = fp_mul(i77,i90)
	i92       = fp_mul(i64,i91)
	i93       = fp_mul(i84,i92)
	i94       = fp_mul(i43,i93)
	i95       = fp_mul(i70,i94)
	i96       = fp_mul(i83,i95)
	i97       = fp_mul(i73,i96)
	i98       = fp_mul(i79,i97)
	i99       = fp_mul(i81,i98)
	i100      = fp_mul(i98,i99)
	i101      = fp_mul(i74,i100)
	i102      = fp_mul(i76,i101)
	i103      = fp_mul(i89,i102)
	i104      = fp_mul(i75,i103)
	i105      = fp_mul(i99,i104)
	i106      = fp_mul(i94,i105)
	i107      = fp_mul(i92,i106)
	i108      = fp_mul(i96,i107)
	i109      = fp_mul(i90,i108)
	i110      = fp_mul(i87,i109)
	i111      = fp_mul(i91,i110)
	i112      = fp_mul(i105,i111)
	i113      = fp_mul(i86,i112)
	i114      = fp_mul(i101,i113)
	i115      = fp_mul(i111,i114)
	i116      = fp_mul(i93,i115)
	i117      = fp_mul(i108,i116)
	i118      = fp_mul(i97,i117)
	i119      = fp_mul(i104,i118)
	i120      = fp_mul(i113,i119)
	i121      = fp_mul(i102,i120)
	i122      = fp_mul(i85,i121)
	i123      = fp_mul(i95,i122)
	i125      = fp_mul(fp_sqr(i123),i107)
	i126      = fp_mul(i121,i125)
	i127      = fp_mul(i120,i126)
	i128      = fp_mul(i114,i127)
	i129      = fp_mul(i115,i128)
	i130      = fp_mul(i112,i129)
	i131      = fp_mul(i117,i130)
	i132      = fp_mul(i100,i131)
	i133      = fp_mul(i110,i132)
	i134      = fp_mul(i119,i133)
	i135      = fp_mul(i109,i134)
	i136      = fp_mul(i123,i135)
	i137      = fp_mul(i106,i136)
	i138      = fp_mul(i116,i137)
	i139      = fp_mul(i118,i138)
	i140      = fp_mul(i103,i139)
	x32       = fp_mul(i122,i140)
	i240      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i136,8589934592),i139),1),4294967296),i134),1),4294967296)
	i306      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i132,i240),1),2147483648),i129),1),4294967296),i126)
	i404      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i306,8589934592),i130),1),8589934592),i137),1),1073741824)
	i475      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i125,i404),1),17179869184),i138),1),17179869184),i140)
	i575      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i475,4294967296),i131),1),34359738368),i135),1),2147483648)
	i647      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i128,i575),1),34359738368),i127),1),17179869184),i133)
	i745      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i647,4294967296),x32),1),4294967296),x32),1),4294967296)
	i812      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i745),1),4294967296),x32),1),4294967296),x32)
	i910      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i812,4294967296),x32),1),4294967296),x32),1),4294967296)
	i977      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i910),1),4294967296),x32),1),4294967296),x32)
	i1075     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i977,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1142     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1075),1),4294967296),x32),1),4294967296),x32)
	i1240     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1142,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1307     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1240),1),4294967296),x32),1),4294967296),x32)
	i1405     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1307,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1472     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1405),1),4294967296),x32),1),4294967296),x32)
	i1570     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1472,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1637     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1570),1),4294967296),x32),1),4294967296),x32)
	i1735     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1637,4294967296),x32),1),4294967296),x32),1),4294967296)
	i1802     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1735),1),4294967296),x32),1),4294967296),x32)
	i1900     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1802,4294967296),x32),1),4294967296),x32),1),4294967296)
	return       fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x32,i1900),1),4294967296),x32),1),1048576),i82)
