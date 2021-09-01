from src.fp import fp_mul, fp_sqr, fp_exp
def quint_4096_3(x):
	_1 = x
	_10       = fp_sqr(_1)
	_100      = fp_sqr(_10)
	_101      = fp_mul(_1,_100)
	_1010     = fp_sqr(_101)
	_1111     = fp_mul(_101,_1010)
	_10001    = fp_mul(_10,_1111)
	_10100    = fp_mul(_101,_1111)
	_11000    = fp_mul(_100,_10100)
	_11100    = fp_mul(_100,_11000)
	_110000   = fp_mul(_10100,_11100)
	_1000001  = fp_mul(_10001,_110000)
	_10000010 = fp_sqr(_1000001)
	_11000011 = fp_mul(_1000001,_10000010)
	_11011011 = fp_mul(_11000,_11000011)
	_11100000 = fp_mul(_101,_11011011)
	i16       = fp_mul(_110000,_11100000)
	i17       = fp_mul(_11100000,i16)
	i18       = fp_mul(_11100,i17)
	i19       = fp_mul(_1000001,i18)
	i20       = fp_mul(i16,i19)
	i21       = fp_mul(_11011011,i20)
	i22       = fp_mul(i19,i21)
	i23       = fp_mul(i17,i22)
	i24       = fp_mul(i18,i23)
	i26       = fp_mul(fp_sqr(i24),i21)
	i27       = fp_mul(i24,i26)
	i28       = fp_mul(i26,i27)
	i29       = fp_mul(i22,i28)
	i30       = fp_mul(i20,i29)
	i31       = fp_mul(i27,i30)
	i32       = fp_mul(i29,i31)
	i33       = fp_mul(i31,i32)
	i34       = fp_mul(i30,i33)
	i35       = fp_mul(i33,i34)
	i36       = fp_mul(i32,i35)
	i37       = fp_mul(i23,i36)
	i39       = fp_mul(fp_sqr(i37),i35)
	i41       = fp_mul(fp_sqr(i39),i36)
	i42       = fp_mul(i28,i41)
	i43       = fp_mul(i39,i42)
	i44       = fp_mul(i37,i43)
	i45       = fp_mul(i34,i44)
	i48       = fp_mul(fp_mul(fp_sqr(i45),i45),i42)
	i50       = fp_mul(fp_sqr(i48),i41)
	i51       = fp_mul(i43,i50)
	i52       = fp_mul(i48,i51)
	i54       = fp_mul(fp_sqr(i52),i45)
	i56       = fp_mul(fp_sqr(i54),i51)
	i58       = fp_mul(fp_sqr(i56),i50)
	i59       = fp_mul(i54,i58)
	i60       = fp_mul(i44,i59)
	i61       = fp_mul(i56,i60)
	i62       = fp_mul(i60,i61)
	i63       = fp_mul(i58,i62)
	i65       = fp_mul(fp_sqr(i63),i52)
	i66       = fp_mul(i61,i65)
	i67       = fp_mul(i59,i66)
	i68       = fp_mul(i63,i67)
	i69       = fp_mul(i62,i68)
	i70       = fp_mul(i68,i69)
	i71       = fp_mul(i67,i70)
	i72       = fp_mul(i69,i71)
	i73       = fp_mul(i65,i72)
	i74       = fp_mul(i72,i73)
	i75       = fp_mul(i66,i74)
	i77       = fp_mul(fp_sqr(i75),i70)
	i78       = fp_mul(i71,i77)
	i80       = fp_mul(fp_sqr(i78),i73)
	i81       = fp_mul(i75,i80)
	i82       = fp_mul(i74,i81)
	i83       = fp_mul(i78,i82)
	i85       = fp_mul(fp_sqr(i83),i77)
	i87       = fp_mul(fp_sqr(i85),i83)
	i88       = fp_mul(i81,i87)
	i89       = fp_mul(i80,i88)
	i90       = fp_mul(i82,i89)
	i91       = fp_mul(i85,i90)
	i92       = fp_mul(i90,i91)
	i93       = fp_mul(i91,i92)
	i94       = fp_mul(i89,i93)
	i95       = fp_mul(i87,i94)
	i96       = fp_mul(i88,i95)
	i97       = fp_mul(i92,i96)
	i98       = fp_mul(i94,i97)
	i99       = fp_mul(i96,i98)
	i100      = fp_mul(i95,i99)
	i101      = fp_mul(i99,i100)
	i102      = fp_mul(i93,i101)
	i103      = fp_mul(i98,i102)
	i104      = fp_mul(i97,i103)
	i105      = fp_mul(i101,i104)
	i106      = fp_mul(i100,i105)
	i107      = fp_mul(i103,i106)
	i109      = fp_mul(fp_sqr(i107),i104)
	i110      = fp_mul(i102,i109)
	i111      = fp_mul(i109,i110)
	i115      = fp_mul(fp_exp(fp_mul(fp_sqr(i111),i111),2),i105)
	i116      = fp_mul(i110,i115)
	i117      = fp_mul(i111,i116)
	i118      = fp_mul(i107,i117)
	i119      = fp_mul(i106,i118)
	i120      = fp_mul(i115,i119)
	i122      = fp_mul(fp_sqr(i120),i118)
	i123      = fp_mul(i116,i122)
	i124      = fp_mul(i122,i123)
	i125      = fp_mul(i119,i124)
	i127      = fp_mul(fp_sqr(i125),i117)
	i128      = fp_mul(i120,i127)
	i129      = fp_mul(i124,i128)
	i130      = fp_mul(i123,i129)
	i131      = fp_mul(i128,i130)
	i132      = fp_mul(i129,i131)
	i133      = fp_mul(i130,i132)
	i134      = fp_mul(i127,i133)
	i135      = fp_mul(i131,i134)
	i136      = fp_mul(i125,i135)
	i137      = fp_mul(i133,i136)
	i138      = fp_mul(i136,i137)
	i140      = fp_mul(fp_sqr(i138),i132)
	i141      = fp_mul(i135,i140)
	i142      = fp_mul(i140,i141)
	i145      = fp_mul(fp_mul(fp_sqr(i142),i142),i134)
	i146      = fp_mul(i137,i145)
	i147      = fp_mul(i142,i146)
	i148      = fp_mul(i146,i147)
	i150      = fp_mul(fp_sqr(i148),i138)
	i151      = fp_mul(i148,i150)
	i152      = fp_mul(i150,i151)
	i153      = fp_mul(i147,i152)
	i154      = fp_mul(i141,i153)
	i155      = fp_mul(i153,i154)
	i157      = fp_mul(fp_sqr(i155),i145)
	i158      = fp_mul(i152,i157)
	i159      = fp_mul(i155,i158)
	i161      = fp_mul(fp_sqr(i159),i151)
	i162      = fp_mul(i158,i161)
	i164      = fp_mul(fp_sqr(i162),i154)
	i165      = fp_mul(i161,i164)
	i166      = fp_mul(i159,i165)
	i167      = fp_mul(i165,i166)
	i168      = fp_mul(i157,i167)
	i170      = fp_mul(fp_sqr(i168),i166)
	i171      = fp_mul(i164,i170)
	i172      = fp_mul(i167,i171)
	i173      = fp_mul(i168,i172)
	i174      = fp_mul(i162,i173)
	i176      = fp_mul(fp_sqr(i174),i170)
	i177      = fp_mul(i171,i176)
	i178      = fp_mul(_1111,i177)
	i180      = fp_mul(fp_sqr(i178),i173)
	i181      = fp_mul(i172,i180)
	i182      = fp_mul(i174,i181)
	i183      = fp_mul(i178,i182)
	i186      = fp_mul(fp_mul(fp_sqr(i183),i183),i181)
	i187      = fp_mul(i177,i186)
	i188      = fp_mul(i180,i187)
	i189      = fp_mul(i183,i188)
	i190      = fp_mul(i176,i189)
	i191      = fp_mul(i187,i190)
	i192      = fp_mul(i182,i191)
	i193      = fp_mul(i189,i192)
	i194      = fp_mul(i191,i193)
	i195      = fp_mul(i193,i194)
	i198      = fp_mul(fp_mul(fp_sqr(i195),i195),i186)
	i199      = fp_mul(i188,i198)
	i200      = fp_mul(i195,i199)
	i201      = fp_mul(i192,i200)
	i202      = fp_mul(i194,i201)
	i203      = fp_mul(i198,i202)
	i204      = fp_mul(i199,i203)
	i205      = fp_mul(i200,i204)
	i206      = fp_mul(i203,i205)
	i207      = fp_mul(i190,i206)
	i208      = fp_mul(i206,i207)
	i209      = fp_mul(i202,i208)
	i210      = fp_mul(i201,i209)
	i211      = fp_mul(i204,i210)
	i212      = fp_mul(i207,i211)
	i213      = fp_mul(i205,i212)
	i214      = fp_mul(i208,i213)
	i215      = fp_mul(i209,i214)
	i216      = fp_mul(i211,i215)
	i217      = fp_mul(i214,i216)
	i218      = fp_mul(i212,i217)
	i219      = fp_mul(i215,i218)
	i220      = fp_mul(i213,i219)
	i221      = fp_mul(i210,i220)
	i222      = fp_mul(i218,i221)
	i223      = fp_mul(i220,i222)
	i224      = fp_mul(i221,i223)
	i225      = fp_mul(i216,i224)
	i227      = fp_mul(fp_sqr(i225),i224)
	i228      = fp_mul(i219,i227)
	i229      = fp_mul(i227,i228)
	i230      = fp_mul(i225,i229)
	i231      = fp_mul(i229,i230)
	i232      = fp_mul(i217,i231)
	i233      = fp_mul(i223,i232)
	i234      = fp_mul(i232,i233)
	i235      = fp_mul(i231,i234)
	i236      = fp_mul(i222,i235)
	i237      = fp_mul(i228,i236)
	i238      = fp_mul(i236,i237)
	i239      = fp_mul(i234,i238)
	i240      = fp_mul(i233,i239)
	i241      = fp_mul(i237,i240)
	i242      = fp_mul(i230,i241)
	i244      = fp_mul(fp_sqr(i242),i238)
	i245      = fp_mul(i242,i244)
	i246      = fp_mul(i239,i245)
	i247      = fp_mul(i240,i246)
	i248      = fp_mul(i235,i247)
	i249      = fp_mul(i246,i248)
	i251      = fp_mul(fp_sqr(i249),i248)
	i254      = fp_mul(fp_exp(i251,4),i244)
	i255      = fp_mul(i249,i254)
	i256      = fp_mul(i247,i255)
	i258      = fp_mul(fp_sqr(i256),i241)
	i259      = fp_mul(i255,i258)
	i260      = fp_mul(i254,i259)
	i261      = fp_mul(i256,i260)
	i262      = fp_mul(i251,i261)
	i263      = fp_mul(i245,i262)
	i264      = fp_mul(i259,i263)
	i266      = fp_mul(fp_sqr(i264),i261)
	i267      = fp_mul(i264,i266)
	i268      = fp_mul(i260,i267)
	i269      = fp_mul(i262,i268)
	i270      = fp_mul(i267,i269)
	i271      = fp_mul(i263,i270)
	i272      = fp_mul(i268,i271)
	i273      = fp_mul(i266,i272)
	i274      = fp_mul(i269,i273)
	x124      = fp_mul(fp_sqr(i274),i271)
	i277      = fp_mul(i258,x124)
	i278      = fp_mul(i273,i277)
	i280      = fp_mul(fp_sqr(i278),i272)
	i281      = fp_mul(i270,i280)
	i282      = fp_mul(i278,i281)
	i283      = fp_mul(i280,i282)
	i284      = fp_mul(i277,i283)
	i285      = fp_mul(i274,i284)
	x128      = fp_mul(i283,i285)
	i675      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i285,170141183460469231731687303715884105728),i282),1),680564733841876926926749214863536422912),i281),1),2722258935367507707706996859454145691648)
	i934      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i284,i675),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1320     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i934,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i1579     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1320),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1965     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1579,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2224     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1965),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i2610     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2224,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2869     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i2610),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i3255     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2869,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i3514     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i3255),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i3900     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i3514,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i4159     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i3900),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	return       fp_mul(fp_exp(i4159,21267647932558653966460912964485513216),x124)