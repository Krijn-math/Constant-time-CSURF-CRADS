from src.fp import fp_mul, fp_sqr, fp_exp
def quint_3072(x):
	_1 = x
	_10       = fp_sqr(_1)
	_100      = fp_sqr(_10)
	_110      = fp_mul(_10,_100)
	_1100     = fp_sqr(_110)
	_1101     = fp_mul(_1,_1100)
	_1111     = fp_mul(_10,_1101)
	_10101    = fp_mul(_110,_1111)
	_11000    = fp_sqr(_1100)
	_100010   = fp_mul(_1101,_10101)
	_101111   = fp_mul(_1101,_100010)
	_111011   = fp_mul(_1100,_101111)
	_1010011  = fp_mul(_11000,_111011)
	_1110101  = fp_mul(_100010,_1010011)
	_11001000 = fp_mul(_1010011,_1110101)
	i15       = fp_mul(_111011,_11001000)
	i16       = fp_mul(_101111,i15)
	i17       = fp_mul(_11001000,i16)
	i19       = fp_mul(fp_sqr(i17),_10101)
	i20       = fp_mul(i16,i19)
	i23       = fp_mul(fp_mul(fp_sqr(i20),i20),i15)
	i24       = fp_mul(i17,i23)
	i25       = fp_mul(_1110101,i24)
	i26       = fp_mul(i23,i25)
	i27       = fp_mul(i25,i26)
	i28       = fp_mul(i26,i27)
	i29       = fp_mul(i27,i28)
	i30       = fp_mul(i24,i29)
	i31       = fp_mul(i19,i30)
	i32       = fp_mul(i28,i31)
	i33       = fp_mul(i31,i32)
	i34       = fp_mul(i30,i33)
	i35       = fp_mul(i33,i34)
	i36       = fp_mul(i32,i35)
	i37       = fp_mul(i20,i36)
	i38       = fp_mul(i34,i37)
	i40       = fp_mul(fp_sqr(i38),i29)
	i43       = fp_mul(fp_mul(fp_sqr(i40),i40),i37)
	i44       = fp_mul(i36,i43)
	i45       = fp_mul(i35,i44)
	i46       = fp_mul(i40,i45)
	i47       = fp_mul(i45,i46)
	i48       = fp_mul(i38,i47)
	i49       = fp_mul(i43,i48)
	i51       = fp_mul(fp_sqr(i49),i44)
	i52       = fp_mul(i47,i51)
	i56       = fp_mul(fp_mul(fp_exp(i52,4),i52),i48)
	i57       = fp_mul(i49,i56)
	i58       = fp_mul(i51,i57)
	i60       = fp_mul(fp_sqr(i58),i56)
	i61       = fp_mul(i52,i60)
	i62       = fp_mul(i46,i61)
	i64       = fp_mul(fp_sqr(i62),i58)
	i66       = fp_mul(fp_sqr(i64),i57)
	i67       = fp_mul(i61,i66)
	i68       = fp_mul(i64,i67)
	i69       = fp_mul(i67,i68)
	i70       = fp_mul(i60,i69)
	i71       = fp_mul(i69,i70)
	i72       = fp_mul(i66,i71)
	i73       = fp_mul(i68,i72)
	i74       = fp_mul(i62,i73)
	i77       = fp_mul(fp_exp(i74,4),i71)
	i78       = fp_mul(i73,i77)
	i80       = fp_mul(fp_sqr(i78),i74)
	i81       = fp_mul(i77,i80)
	i82       = fp_mul(i78,i81)
	i83       = fp_mul(i70,i82)
	i84       = fp_mul(i80,i83)
	i85       = fp_mul(i81,i84)
	i86       = fp_mul(i72,i85)
	i89       = fp_mul(fp_mul(fp_sqr(i86),i86),i82)
	i90       = fp_mul(i85,i89)
	i91       = fp_mul(i84,i90)
	i92       = fp_mul(i90,i91)
	i93       = fp_mul(i86,i92)
	i94       = fp_mul(i83,i93)
	i97       = fp_mul(fp_mul(fp_sqr(i94),i94),i91)
	i98       = fp_mul(i93,i97)
	i99       = fp_mul(i89,i98)
	i100      = fp_mul(i97,i99)
	i101      = fp_mul(i99,i100)
	i102      = fp_mul(i98,i101)
	i103      = fp_mul(i92,i102)
	i104      = fp_mul(i94,i103)
	i105      = fp_mul(i103,i104)
	i106      = fp_mul(i101,i105)
	i107      = fp_mul(i100,i106)
	i108      = fp_mul(i102,i107)
	i109      = fp_mul(i104,i108)
	i110      = fp_mul(i105,i109)
	i111      = fp_mul(i106,i110)
	i112      = fp_mul(i107,i111)
	i114      = fp_mul(fp_sqr(i112),i110)
	i115      = fp_mul(i112,i114)
	i116      = fp_mul(i108,i115)
	i118      = fp_mul(fp_sqr(i116),i115)
	i119      = fp_mul(i114,i118)
	i120      = fp_mul(i116,i119)
	i121      = fp_mul(i118,i120)
	i124      = fp_mul(fp_mul(fp_sqr(i121),i121),i111)
	i125      = fp_mul(i121,i124)
	i127      = fp_mul(fp_sqr(i125),i119)
	i128      = fp_mul(i109,i127)
	i129      = fp_mul(i120,i128)
	i130      = fp_mul(i125,i129)
	i131      = fp_mul(i127,i130)
	i133      = fp_mul(fp_sqr(i131),i128)
	i134      = fp_mul(i131,i133)
	i135      = fp_mul(i130,i134)
	i136      = fp_mul(i124,i135)
	i137      = fp_mul(i129,i136)
	i139      = fp_mul(fp_sqr(i137),i134)
	i140      = fp_mul(i133,i139)
	i142      = fp_mul(fp_sqr(i140),i136)
	i143      = fp_mul(i140,i142)
	i144      = fp_mul(i135,i143)
	i145      = fp_mul(i139,i144)
	i146      = fp_mul(i137,i145)
	i148      = fp_mul(fp_sqr(i146),i143)
	i149      = fp_mul(i145,i148)
	i151      = fp_mul(fp_sqr(i149),i146)
	i152      = fp_mul(i144,i151)
	i153      = fp_mul(i142,i152)
	i154      = fp_mul(i149,i153)
	i155      = fp_mul(i148,i154)
	i159      = fp_mul(fp_mul(fp_exp(i155,4),i155),i151)
	i160      = fp_mul(i155,i159)
	i161      = fp_mul(i152,i160)
	i163      = fp_mul(fp_sqr(i161),i159)
	i164      = fp_mul(i160,i163)
	i165      = fp_mul(i154,i164)
	i167      = fp_mul(fp_sqr(i165),i163)
	i168      = fp_mul(i164,i167)
	i169      = fp_mul(i161,i168)
	i171      = fp_mul(fp_sqr(i169),i168)
	i172      = fp_mul(i153,i171)
	i173      = fp_mul(i167,i172)
	i174      = fp_mul(i172,i173)
	i175      = fp_mul(i173,i174)
	i176      = fp_mul(i165,i175)
	i177      = fp_mul(i171,i176)
	i179      = fp_mul(fp_sqr(i177),i169)
	i180      = fp_mul(i174,i179)
	i181      = fp_mul(i175,i180)
	i182      = fp_mul(i179,i181)
	i183      = fp_mul(i177,i182)
	i184      = fp_mul(i182,i183)
	i185      = fp_mul(i176,i184)
	i186      = fp_mul(i183,i185)
	i187      = fp_mul(i181,i186)
	i188      = fp_mul(i184,i187)
	i189      = fp_mul(i180,i188)
	i190      = fp_mul(i185,i189)
	i191      = fp_mul(i187,i190)
	i192      = fp_mul(i189,i191)
	i193      = fp_mul(i186,i192)
	i195      = fp_mul(fp_sqr(i193),i192)
	i196      = fp_mul(i191,i195)
	i199      = fp_mul(fp_exp(i196,4),i195)
	i200      = fp_mul(i190,i199)
	i201      = fp_mul(i196,i200)
	i202      = fp_mul(i193,i201)
	i203      = fp_mul(i200,i202)
	i204      = fp_mul(i188,i203)
	i205      = fp_mul(i199,i204)
	i206      = fp_mul(i201,i205)
	i207      = fp_mul(i204,i206)
	i208      = fp_mul(i202,i207)
	i209      = fp_mul(i206,i208)
	i210      = fp_mul(i203,i209)
	i211      = fp_mul(i208,i210)
	i212      = fp_mul(i209,i211)
	i213      = fp_mul(i210,i212)
	i214      = fp_mul(i207,i213)
	i215      = fp_mul(i205,i214)
	i216      = fp_mul(i211,i215)
	i217      = fp_mul(i213,i216)
	i218      = fp_mul(i212,i217)
	i220      = fp_mul(fp_sqr(i218),i214)
	i221      = fp_mul(i218,i220)
	i222      = fp_mul(_1111,i221)
	i223      = fp_mul(i220,i222)
	i224      = fp_mul(i217,i223)
	i225      = fp_mul(i221,i224)
	i226      = fp_mul(i215,i225)
	i227      = fp_mul(i223,i226)
	i228      = fp_mul(i225,i227)
	i229      = fp_mul(i226,i228)
	i230      = fp_mul(i224,i229)
	i231      = fp_mul(i222,i230)
	i233      = fp_mul(fp_sqr(i231),i216)
	i234      = fp_mul(i230,i233)
	i235      = fp_mul(i233,i234)
	i236      = fp_mul(i229,i235)
	i239      = fp_mul(fp_mul(fp_sqr(i236),i236),i228)
	i240      = fp_mul(i235,i239)
	i241      = fp_mul(i236,i240)
	i242      = fp_mul(i227,i241)
	i243      = fp_mul(i231,i242)
	i246      = fp_mul(fp_mul(fp_sqr(i243),i243),i234)
	i247      = fp_mul(i241,i246)
	i248      = fp_mul(i243,i247)
	i249      = fp_mul(i239,i248)
	i250      = fp_mul(i240,i249)
	i251      = fp_mul(i249,i250)
	i252      = fp_mul(i247,i251)
	i254      = fp_mul(fp_sqr(i252),i246)
	i255      = fp_mul(i250,i254)
	i256      = fp_mul(i248,i255)
	i259      = fp_mul(fp_mul(fp_sqr(i256),i256),i254)
	i260      = fp_mul(i242,i259)
	i261      = fp_mul(i251,i260)
	i262      = fp_mul(i260,i261)
	i263      = fp_mul(i255,i262)
	i264      = fp_mul(i259,i263)
	i265      = fp_mul(i262,i264)
	i266      = fp_mul(i261,i265)
	i267      = fp_mul(i256,i266)
	i268      = fp_mul(i252,i267)
	i269      = fp_mul(i263,i268)
	i270      = fp_mul(i268,i269)
	i271      = fp_mul(i265,i270)
	i272      = fp_mul(i270,i271)
	i273      = fp_mul(i264,i272)
	i274      = fp_mul(i269,i273)
	x124      = fp_mul(i266,i274)
	i276      = fp_mul(i271,x124)
	i278      = fp_mul(fp_sqr(i276),i267)
	i279      = fp_mul(i272,i278)
	i280      = fp_mul(i274,i279)
	i281      = fp_mul(i278,i280)
	i283      = fp_mul(fp_mul(i273,i281),i279)
	x128      = fp_mul(i280,i283)
	i542      = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i264,x128),1),340282366920938463463374607431768211456),i281),1),170141183460469231731687303715884105728),i276)
	i931      = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i542,2722258935367507707706996859454145691648),i283),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i1190     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i931),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1576     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1190,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i1835     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1576),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i2221     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1835,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2480     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i2221),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i2866     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2480,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i3125     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i2866),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	return       fp_mul(fp_exp(i3125,21267647932558653966460912964485513216),x124)
