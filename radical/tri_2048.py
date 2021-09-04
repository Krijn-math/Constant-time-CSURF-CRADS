from src.fp import fp_mul, fp_sqr, fp_exp
def tri_2048(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_110      = fp_sqr(_11)
	_1100     = fp_sqr(_110)
	_1101     = fp_mul(_1,_1100)
	_1110     = fp_mul(_1,_1101)
	_10000    = fp_mul(_10,_1110)
	_11110    = fp_mul(_1110,_10000)
	_100001   = fp_mul(_11,_11110)
	_100100   = fp_mul(_11,_100001)
	_100101   = fp_mul(_1,_100100)
	_1000110  = fp_mul(_100001,_100101)
	_1010110  = fp_mul(_10000,_1000110)
	_1011001  = fp_mul(_11,_1010110)
	_1011111  = fp_mul(_110,_1011001)
	_10000100 = fp_mul(_100101,_1011111)
	_10100010 = fp_mul(_11110,_10000100)
	_11000110 = fp_mul(_100100,_10100010)
	_11010010 = fp_mul(_1100,_11000110)
	i20       = fp_mul(_11000110,_11010010)
	i21       = fp_mul(_10000100,i20)
	i22       = fp_mul(_1010110,i21)
	i23       = fp_mul(_1000110,i22)
	i24       = fp_mul(i21,i23)
	i25       = fp_mul(i23,i24)
	i26       = fp_mul(_1011001,i25)
	i27       = fp_mul(i20,i26)
	i28       = fp_mul(_11010010,i27)
	i29       = fp_mul(_1101,i28)
	i30       = fp_mul(_10100010,i29)
	i31       = fp_mul(i22,i30)
	i32       = fp_mul(_1011001,i31)
	i33       = fp_mul(i30,i32)
	i34       = fp_mul(_1011111,i33)
	i35       = fp_mul(i26,i34)
	i36       = fp_mul(i31,i35)
	i37       = fp_mul(i24,i36)
	i38       = fp_mul(i34,i37)
	i39       = fp_mul(i27,i38)
	i40       = fp_mul(i33,i39)
	i41       = fp_mul(i37,i40)
	i42       = fp_mul(i29,i41)
	i43       = fp_mul(i36,i42)
	i44       = fp_mul(i39,i43)
	i45       = fp_mul(i32,i44)
	i46       = fp_mul(i25,i45)
	i47       = fp_mul(i38,i46)
	i48       = fp_mul(i44,i47)
	i49       = fp_mul(i40,i48)
	i50       = fp_mul(i35,i49)
	i51       = fp_mul(i46,i50)
	i52       = fp_mul(i50,i51)
	i53       = fp_mul(i43,i52)
	i54       = fp_mul(i41,i53)
	i55       = fp_mul(i45,i54)
	i56       = fp_mul(i49,i55)
	i57       = fp_mul(i47,i56)
	i58       = fp_mul(i28,i57)
	i59       = fp_mul(i55,i58)
	i60       = fp_mul(i57,i59)
	i61       = fp_mul(i53,i60)
	i62       = fp_mul(i56,i61)
	i63       = fp_mul(i54,i62)
	i64       = fp_mul(i59,i63)
	i65       = fp_mul(i51,i64)
	i66       = fp_mul(i58,i65)
	i67       = fp_mul(i65,i66)
	i68       = fp_mul(i60,i67)
	i69       = fp_mul(i61,i68)
	i70       = fp_mul(i63,i69)
	i71       = fp_mul(i52,i70)
	i72       = fp_mul(i48,i71)
	i73       = fp_mul(i64,i72)
	i74       = fp_mul(i69,i73)
	i75       = fp_mul(i62,i74)
	i76       = fp_mul(i66,i75)
	i77       = fp_mul(i73,i76)
	i78       = fp_mul(i70,i77)
	i79       = fp_mul(i67,i78)
	i80       = fp_mul(i74,i79)
	i81       = fp_mul(i71,i80)
	i82       = fp_mul(i72,i81)
	i83       = fp_mul(i75,i82)
	i84       = fp_mul(i77,i83)
	i85       = fp_mul(i68,i84)
	i86       = fp_mul(i80,i85)
	i87       = fp_mul(i82,i86)
	i88       = fp_mul(i81,i87)
	i89       = fp_mul(i78,i88)
	i90       = fp_mul(i42,i89)
	i91       = fp_mul(i76,i90)
	i92       = fp_mul(i83,i91)
	i93       = fp_mul(i85,i92)
	i94       = fp_mul(i86,i93)
	i95       = fp_mul(i84,i94)
	i96       = fp_mul(i94,i95)
	i97       = fp_mul(i89,i96)
	i98       = fp_mul(i79,i97)
	i99       = fp_mul(i88,i98)
	i101      = fp_mul(fp_sqr(i99),i91)
	i102      = fp_mul(i99,i101)
	i103      = fp_mul(i93,i102)
	i104      = fp_mul(i98,i103)
	i105      = fp_mul(i90,i104)
	i106      = fp_mul(i102,i105)
	i107      = fp_mul(i97,i106)
	i108      = fp_mul(i101,i107)
	i109      = fp_mul(i95,i108)
	i110      = fp_mul(i96,i109)
	i111      = fp_mul(i105,i110)
	i112      = fp_mul(i103,i111)
	i113      = fp_mul(i108,i112)
	i114      = fp_mul(i107,i113)
	i115      = fp_mul(i109,i114)
	i116      = fp_mul(i87,i115)
	i118      = fp_mul(fp_sqr(i116),i92)
	i119      = fp_mul(i104,i118)
	i120      = fp_mul(i113,i119)
	i121      = fp_mul(i110,i120)
	i122      = fp_mul(i115,i121)
	i123      = fp_mul(i114,i122)
	i124      = fp_mul(i106,i123)
	i125      = fp_mul(i112,i124)
	i126      = fp_mul(i121,i125)
	i127      = fp_mul(i118,i126)
	i128      = fp_mul(i120,i127)
	i129      = fp_mul(i122,i128)
	i131      = fp_mul(fp_sqr(i129),i119)
	i132      = fp_mul(i111,i131)
	i133      = fp_mul(i126,i132)
	i134      = fp_mul(i123,i133)
	i135      = fp_mul(i125,i134)
	i136      = fp_mul(i132,i135)
	i137      = fp_mul(i128,i136)
	i138      = fp_mul(i129,i137)
	i139      = fp_mul(i131,i138)
	i140      = fp_mul(i124,i139)
	i141      = fp_mul(i136,i140)
	i142      = fp_mul(i127,i141)
	i143      = fp_mul(i137,i142)
	i144      = fp_mul(i135,i143)
	i145      = fp_mul(i133,i144)
	i146      = fp_mul(i116,i145)
	i147      = fp_mul(i134,i146)
	i148      = fp_mul(i138,i147)
	i149      = fp_mul(i140,i148)
	i150      = fp_mul(i139,i149)
	i151      = fp_mul(i143,i150)
	i152      = fp_mul(i142,i151)
	i153      = fp_mul(i146,i152)
	i154      = fp_mul(i147,i153)
	i155      = fp_mul(i141,i154)
	i157      = fp_mul(fp_sqr(i155),i145)
	i158      = fp_mul(i155,i157)
	i159      = fp_mul(i144,i158)
	i160      = fp_mul(i157,i159)
	i161      = fp_mul(i152,i160)
	i162      = fp_mul(i159,i161)
	i163      = fp_mul(i154,i162)
	i164      = fp_mul(i148,i163)
	i166      = fp_mul(fp_sqr(i164),i149)
	i167      = fp_mul(i150,i166)
	i168      = fp_mul(i162,i167)
	i169      = fp_mul(i151,i168)
	i170      = fp_mul(i153,i169)
	i171      = fp_mul(i160,i170)
	i172      = fp_mul(i161,i171)
	i173      = fp_mul(i158,i172)
	i174      = fp_mul(i168,i173)
	i175      = fp_mul(i164,i174)
	i176      = fp_mul(i169,i175)
	x57       = fp_mul(i163,i176)
	i178      = fp_mul(i166,x57)
	i179      = fp_mul(i173,i178)
	i180      = fp_mul(i167,i179)
	i181      = fp_mul(i172,i180)
	i182      = fp_mul(i180,i181)
	i183      = fp_mul(i176,i182)
	i184      = fp_mul(i178,i183)
	i185      = fp_mul(i171,i184)
	i186      = fp_mul(i170,i185)
	i187      = fp_mul(i175,i186)
	i188      = fp_mul(i181,i187)
	i189      = fp_mul(i174,i188)
	i190      = fp_mul(i179,i189)
	i191      = fp_mul(i188,i190)
	i192      = fp_mul(i186,i191)
	i193      = fp_mul(i185,i192)
	i194      = fp_mul(i182,i193)
	i195      = fp_mul(i183,i194)
	i196      = fp_mul(i192,i195)
	i197      = fp_mul(i190,i196)
	i198      = fp_mul(i193,i197)
	i199      = fp_mul(i189,i198)
	x64       = fp_mul(i187,i199)
	i396      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i196,18446744073709551616),i198),1),9223372036854775808),i184),1),147573952589676412928)
	i530      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i191,i396),1),295147905179352825856),i197),1),9223372036854775808),i195)
	i725      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i530,18446744073709551616),i194),1),36893488147419103232),i199),1),18446744073709551616)
	i856      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i725),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i1050     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i856,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i1181     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i1050),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i1375     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1181,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i1506     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i1375),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i1700     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1506,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i1831     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i1700),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i2025     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1831,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i2156     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i2025),1),18446744073709551616),x64),1),18446744073709551616),x64)
	return       fp_mul(fp_exp(i2156,144115188075855872),x57)
