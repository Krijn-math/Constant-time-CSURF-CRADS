from src.fp import fp_mul, fp_sqr, fp_exp
def quint_2048(x):
	_1 = x
	_10       = fp_sqr(_1)
	_11       = fp_mul(_1,_10)
	_110      = fp_sqr(_11)
	_1000     = fp_mul(_10,_110)
	_1001     = fp_mul(_1,_1000)
	_1010     = fp_mul(_1,_1001)
	_10001    = fp_mul(_1000,_1001)
	_10010    = fp_mul(_1,_10001)
	_100100   = fp_sqr(_10010)
	_110110   = fp_mul(_10010,_100100)
	_1000111  = fp_mul(_10001,_110110)
	_1001111  = fp_mul(_1000,_1000111)
	_1011010  = fp_mul(_100100,_110110)
	_1100100  = fp_mul(_1010,_1011010)
	_1100101  = fp_mul(_1,_1100100)
	_1110110  = fp_mul(_10001,_1100101)
	_1111111  = fp_mul(_1001,_1110110)
	_10001110 = fp_sqr(_1000111)
	_10011010 = fp_mul(_100100,_1110110)
	_11100001 = fp_mul(_1000111,_10011010)
	i21       = fp_mul(_1001111,_11100001)
	i22       = fp_mul(_10010,i21)
	i23       = fp_mul(_1110110,i22)
	i24       = fp_mul(_1100100,i23)
	i25       = fp_mul(i23,i24)
	i26       = fp_mul(_1100101,i25)
	i28       = fp_mul(fp_sqr(i26),_11100001)
	i29       = fp_mul(_1111111,i28)
	i30       = fp_mul(i21,i29)
	i31       = fp_mul(i24,i30)
	i32       = fp_mul(i25,i31)
	i33       = fp_mul(i26,i32)
	i34       = fp_mul(_10001110,i33)
	i35       = fp_mul(i28,i34)
	i36       = fp_mul(_10011010,i35)
	i37       = fp_mul(i32,i36)
	i38       = fp_mul(i29,i37)
	i39       = fp_mul(i31,i38)
	i40       = fp_mul(i22,i39)
	i41       = fp_mul(i30,i40)
	i42       = fp_mul(_1011010,i41)
	i43       = fp_mul(i37,i42)
	i44       = fp_mul(i42,i43)
	i45       = fp_mul(i43,i44)
	i46       = fp_mul(i34,i45)
	i47       = fp_mul(i35,i46)
	i48       = fp_mul(i33,i47)
	i49       = fp_mul(i47,i48)
	i50       = fp_mul(i45,i49)
	i51       = fp_mul(i39,i50)
	i52       = fp_mul(i46,i51)
	i53       = fp_mul(i44,i52)
	i54       = fp_mul(i49,i53)
	i55       = fp_mul(i48,i54)
	i56       = fp_mul(i40,i55)
	i57       = fp_mul(i36,i56)
	i58       = fp_mul(i56,i57)
	i59       = fp_mul(i51,i58)
	i60       = fp_mul(i52,i59)
	i61       = fp_mul(i38,i60)
	i62       = fp_mul(i57,i61)
	i63       = fp_mul(i54,i62)
	i64       = fp_mul(i41,i63)
	i65       = fp_mul(i59,i64)
	i66       = fp_mul(i50,i65)
	i67       = fp_mul(i63,i66)
	i68       = fp_mul(i64,i67)
	i69       = fp_mul(i53,i68)
	i70       = fp_mul(i65,i69)
	i71       = fp_mul(i66,i70)
	i72       = fp_mul(i55,i71)
	i73       = fp_mul(i71,i72)
	i74       = fp_mul(i58,i73)
	i75       = fp_mul(i67,i74)
	i76       = fp_mul(i61,i75)
	i77       = fp_mul(i60,i76)
	i78       = fp_mul(i68,i77)
	i79       = fp_mul(i62,i78)
	i80       = fp_mul(i69,i79)
	i81       = fp_mul(i78,i80)
	i82       = fp_mul(i76,i81)
	i83       = fp_mul(i77,i82)
	i84       = fp_mul(i74,i83)
	i85       = fp_mul(i79,i84)
	i86       = fp_mul(i75,i85)
	i87       = fp_mul(i81,i86)
	i88       = fp_mul(i72,i87)
	i89       = fp_mul(i85,i88)
	i90       = fp_mul(i84,i89)
	i91       = fp_mul(i82,i90)
	i92       = fp_mul(i70,i91)
	i93       = fp_mul(i83,i92)
	i94       = fp_mul(i73,i93)
	i95       = fp_mul(i91,i94)
	i96       = fp_mul(i87,i95)
	i97       = fp_mul(i86,i96)
	i98       = fp_mul(i94,i97)
	i99       = fp_mul(i90,i98)
	i100      = fp_mul(i93,i99)
	i101      = fp_mul(i80,i100)
	i103      = fp_mul(fp_sqr(i101),i88)
	i104      = fp_mul(i89,i103)
	i105      = fp_mul(i96,i104)
	i106      = fp_mul(i97,i105)
	i107      = fp_mul(i101,i106)
	i108      = fp_mul(i99,i107)
	i109      = fp_mul(i98,i108)
	i110      = fp_mul(i95,i109)
	i111      = fp_mul(i92,i110)
	i112      = fp_mul(i103,i111)
	i113      = fp_mul(i105,i112)
	i114      = fp_mul(i111,i113)
	i115      = fp_mul(i107,i114)
	i116      = fp_mul(i109,i115)
	i117      = fp_mul(i110,i116)
	i118      = fp_mul(i108,i117)
	i120      = fp_mul(fp_sqr(i118),i100)
	i121      = fp_mul(i117,i120)
	i122      = fp_mul(i118,i121)
	i123      = fp_mul(i115,i122)
	i124      = fp_mul(i112,i123)
	i125      = fp_mul(i116,i124)
	i126      = fp_mul(i120,i125)
	i127      = fp_mul(i125,i126)
	i128      = fp_mul(i104,i127)
	i129      = fp_mul(i123,i128)
	i130      = fp_mul(i121,i129)
	i131      = fp_mul(i106,i130)
	i132      = fp_mul(i113,i131)
	i134      = fp_mul(fp_sqr(i132),i122)
	i135      = fp_mul(i131,i134)
	i136      = fp_mul(i114,i135)
	i137      = fp_mul(i130,i136)
	i138      = fp_mul(i126,i137)
	i139      = fp_mul(i128,i138)
	i140      = fp_mul(i132,i139)
	i141      = fp_mul(i134,i140)
	i142      = fp_mul(i129,i141)
	i144      = fp_mul(fp_sqr(i142),i139)
	i145      = fp_mul(i124,i144)
	i146      = fp_mul(i136,i145)
	i147      = fp_mul(i135,i146)
	i148      = fp_mul(i137,i147)
	i149      = fp_mul(i147,i148)
	i150      = fp_mul(i127,i149)
	i151      = fp_mul(i140,i150)
	i152      = fp_mul(i142,i151)
	i153      = fp_mul(i138,i152)
	i154      = fp_mul(i146,i153)
	i156      = fp_mul(fp_sqr(i154),i149)
	i157      = fp_mul(i152,i156)
	i158      = fp_mul(i148,i157)
	i159      = fp_mul(i141,i158)
	i160      = fp_mul(i150,i159)
	i161      = fp_mul(i154,i160)
	i162      = fp_mul(i151,i161)
	i163      = fp_mul(i144,i162)
	i164      = fp_mul(i153,i163)
	i165      = fp_mul(i158,i164)
	i166      = fp_mul(i164,i165)
	i167      = fp_mul(i160,i166)
	i168      = fp_mul(i163,i167)
	i169      = fp_mul(i157,i168)
	i170      = fp_mul(i166,i169)
	i171      = fp_mul(i156,i170)
	i172      = fp_mul(i145,i171)
	i173      = fp_mul(i165,i172)
	i174      = fp_mul(i162,i173)
	i175      = fp_mul(i170,i174)
	i176      = fp_mul(i173,i175)
	i177      = fp_mul(i169,i176)
	i178      = fp_mul(i172,i177)
	i179      = fp_mul(i161,i178)
	i180      = fp_mul(i168,i179)
	x58       = fp_mul(i159,i180)
	i182      = fp_mul(i167,x58)
	i183      = fp_mul(i171,i182)
	i184      = fp_mul(i177,i183)
	i185      = fp_mul(i178,i184)
	i186      = fp_mul(i175,i185)
	i187      = fp_mul(i184,i186)
	i188      = fp_mul(i180,i187)
	i189      = fp_mul(i183,i188)
	i190      = fp_mul(i176,i189)
	i191      = fp_mul(i179,i190)
	i192      = fp_mul(i174,i191)
	i193      = fp_mul(i185,i192)
	i194      = fp_mul(i187,i193)
	i195      = fp_mul(i189,i194)
	i196      = fp_mul(i192,i195)
	i197      = fp_mul(i182,i196)
	i198      = fp_mul(i188,i197)
	i199      = fp_mul(i190,i198)
	i200      = fp_mul(i193,i199)
	x64       = fp_mul(i186,i200)
	i396      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i200,9223372036854775808),i195),1),36893488147419103232),i199),1),36893488147419103232)
	i530      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i197,i396),1),9223372036854775808),i194),1),295147905179352825856),i198)
	i724      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i530,4611686018427387904),i191),1),73786976294838206464),i196),1),18446744073709551616)
	i855      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i724),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i1049     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i855,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i1180     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i1049),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i1374     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1180,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i1505     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i1374),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i1699     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1505,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i1830     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i1699),1),18446744073709551616),x64),1),18446744073709551616),x64)
	i2024     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1830,18446744073709551616),x64),1),18446744073709551616),x64),1),18446744073709551616)
	i2155     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x64,i2024),1),18446744073709551616),x64),1),18446744073709551616),x64)
	return       fp_mul(fp_exp(i2155,288230376151711744),x58)
