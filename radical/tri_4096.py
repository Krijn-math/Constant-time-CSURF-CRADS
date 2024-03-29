from src.fp import fp_mul, fp_sqr, fp_exp
def tri_4096(x):
	_1 = x
	_10      = fp_sqr(_1)
	_100     = fp_sqr(_10)
	_101     = fp_mul(_1,_100)
	_1010    = fp_sqr(_101)
	_1011    = fp_mul(_1,_1010)
	_1101    = fp_mul(_10,_1011)
	_11010   = fp_sqr(_1101)
	_11111   = fp_mul(_101,_11010)
	_100011  = fp_mul(_100,_11111)
	_100101  = fp_mul(_10,_100011)
	_110000  = fp_mul(_1011,_100101)
	_1100000 = fp_sqr(_110000)
	_1100001 = fp_mul(_1,_1100000)
	i17      = fp_mul(fp_mul(fp_exp(_1100001,4),_1100001),_1101)
	i18      = fp_mul(_110000,i17)
	i19      = fp_mul(_1100001,i18)
	i20      = fp_mul(_100101,i19)
	i21      = fp_mul(_100011,i20)
	i22      = fp_mul(i20,i21)
	i23      = fp_mul(i17,i22)
	i24      = fp_mul(i22,i23)
	i25      = fp_mul(i19,i24)
	i26      = fp_mul(i18,i25)
	i27      = fp_mul(i23,i26)
	i28      = fp_mul(i24,i27)
	i29      = fp_mul(i26,i28)
	i31      = fp_mul(fp_sqr(i29),i25)
	i32      = fp_mul(i28,i31)
	i33      = fp_mul(i27,i32)
	i35      = fp_mul(fp_sqr(i33),i29)
	i36      = fp_mul(i21,i35)
	i37      = fp_mul(i32,i36)
	i38      = fp_mul(i33,i37)
	i39      = fp_mul(i31,i38)
	i40      = fp_mul(i35,i39)
	i43      = fp_mul(fp_mul(fp_sqr(i40),i40),i38)
	i44      = fp_mul(i36,i43)
	i45      = fp_mul(i39,i44)
	i48      = fp_mul(fp_exp(i45,4),i44)
	i49      = fp_mul(i40,i48)
	i50      = fp_mul(i45,i49)
	i51      = fp_mul(i37,i50)
	i53      = fp_mul(fp_sqr(i51),i43)
	i54      = fp_mul(i48,i53)
	i55      = fp_mul(i53,i54)
	i56      = fp_mul(i54,i55)
	i58      = fp_mul(fp_sqr(i56),i50)
	i59      = fp_mul(i56,i58)
	i60      = fp_mul(i51,i59)
	i61      = fp_mul(i49,i60)
	i62      = fp_mul(i59,i61)
	i64      = fp_mul(fp_sqr(i62),i61)
	i65      = fp_mul(i60,i64)
	i66      = fp_mul(i55,i65)
	i67      = fp_mul(i62,i66)
	i70      = fp_mul(fp_exp(i67,4),i65)
	i71      = fp_mul(i58,i70)
	i72      = fp_mul(i66,i71)
	i73      = fp_mul(i64,i72)
	i74      = fp_mul(i67,i73)
	i76      = fp_mul(fp_sqr(i74),i71)
	i77      = fp_mul(i73,i76)
	i79      = fp_mul(fp_sqr(i77),i70)
	i80      = fp_mul(i77,i79)
	i81      = fp_mul(i74,i80)
	i83      = fp_mul(fp_sqr(i81),i79)
	i84      = fp_mul(i80,i83)
	i85      = fp_mul(i72,i84)
	i86      = fp_mul(i83,i85)
	i87      = fp_mul(i76,i86)
	i88      = fp_mul(i81,i87)
	i91      = fp_mul(fp_mul(fp_sqr(i88),i88),i86)
	i94      = fp_mul(fp_mul(fp_sqr(i91),i91),i88)
	i95      = fp_mul(i85,i94)
	i96      = fp_mul(i87,i95)
	i97      = fp_mul(i94,i96)
	i98      = fp_mul(i95,i97)
	i99      = fp_mul(i84,i98)
	i100     = fp_mul(i98,i99)
	i101     = fp_mul(i97,i100)
	i102     = fp_mul(i100,i101)
	i103     = fp_mul(i91,i102)
	i104     = fp_mul(i99,i103)
	i105     = fp_mul(i96,i104)
	i106     = fp_mul(i102,i105)
	i107     = fp_mul(i103,i106)
	i108     = fp_mul(i101,i107)
	i109     = fp_mul(i104,i108)
	i110     = fp_mul(i107,i109)
	i111     = fp_mul(i105,i110)
	i112     = fp_mul(i108,i111)
	i113     = fp_mul(i109,i112)
	i114     = fp_mul(i106,i113)
	i115     = fp_mul(i110,i114)
	i116     = fp_mul(i111,i115)
	i117     = fp_mul(i113,i116)
	i119     = fp_mul(fp_sqr(i117),i114)
	i120     = fp_mul(i112,i119)
	i122     = fp_mul(fp_sqr(i120),i117)
	i123     = fp_mul(i116,i122)
	i124     = fp_mul(i119,i123)
	i125     = fp_mul(i115,i124)
	i128     = fp_mul(fp_mul(fp_sqr(i125),i125),i120)
	i129     = fp_mul(i123,i128)
	i131     = fp_mul(fp_sqr(i129),i122)
	i132     = fp_mul(i124,i131)
	i135     = fp_mul(fp_mul(fp_sqr(i132),i132),i125)
	i136     = fp_mul(i132,i135)
	i137     = fp_mul(i131,i136)
	i138     = fp_mul(i135,i137)
	i139     = fp_mul(i137,i138)
	i141     = fp_mul(fp_sqr(i139),i136)
	i142     = fp_mul(i128,i141)
	i143     = fp_mul(i141,i142)
	i144     = fp_mul(i138,i143)
	i145     = fp_mul(i142,i144)
	i146     = fp_mul(i129,i145)
	i148     = fp_mul(fp_sqr(i146),i139)
	i149     = fp_mul(i146,i148)
	i150     = fp_mul(i143,i149)
	i151     = fp_mul(i144,i150)
	i154     = fp_mul(fp_mul(fp_sqr(i151),i151),i145)
	i155     = fp_mul(i148,i154)
	i156     = fp_mul(i154,i155)
	i157     = fp_mul(i151,i156)
	i158     = fp_mul(i150,i157)
	i160     = fp_mul(fp_sqr(i158),i155)
	i161     = fp_mul(i157,i160)
	i162     = fp_mul(i160,i161)
	i163     = fp_mul(i149,i162)
	i164     = fp_mul(i158,i163)
	i165     = fp_mul(i163,i164)
	i166     = fp_mul(i156,i165)
	i168     = fp_mul(fp_sqr(i166),i162)
	i169     = fp_mul(i161,i168)
	i170     = fp_mul(i164,i169)
	i174     = fp_mul(fp_exp(fp_mul(fp_sqr(i170),i170),2),i165)
	i175     = fp_mul(_11111,i174)
	i176     = fp_mul(i168,i175)
	i178     = fp_mul(fp_sqr(i176),i166)
	i181     = fp_mul(fp_mul(fp_sqr(i178),i178),i169)
	i182     = fp_mul(i178,i181)
	i183     = fp_mul(i176,i182)
	i184     = fp_mul(i182,i183)
	i187     = fp_mul(fp_mul(fp_sqr(i184),i184),i183)
	i188     = fp_mul(i184,i187)
	i189     = fp_mul(i175,i188)
	i190     = fp_mul(i181,i189)
	i191     = fp_mul(i188,i190)
	i193     = fp_mul(fp_sqr(i191),i170)
	i194     = fp_mul(i174,i193)
	i195     = fp_mul(i187,i194)
	i196     = fp_mul(i190,i195)
	i197     = fp_mul(i195,i196)
	i198     = fp_mul(i191,i197)
	i199     = fp_mul(i194,i198)
	i200     = fp_mul(i197,i199)
	i201     = fp_mul(i193,i200)
	i202     = fp_mul(i200,i201)
	i203     = fp_mul(i198,i202)
	i204     = fp_mul(i189,i203)
	i205     = fp_mul(i196,i204)
	i206     = fp_mul(i201,i205)
	i207     = fp_mul(i203,i206)
	i208     = fp_mul(i202,i207)
	i210     = fp_mul(fp_sqr(i208),i205)
	i211     = fp_mul(i204,i210)
	i212     = fp_mul(i199,i211)
	i213     = fp_mul(i207,i212)
	i214     = fp_mul(i208,i213)
	i215     = fp_mul(i211,i214)
	i217     = fp_mul(fp_sqr(i215),i210)
	i218     = fp_mul(i206,i217)
	i219     = fp_mul(i212,i218)
	i220     = fp_mul(i215,i219)
	i221     = fp_mul(i213,i220)
	i222     = fp_mul(i219,i221)
	i223     = fp_mul(i214,i222)
	i224     = fp_mul(i220,i223)
	i225     = fp_mul(i217,i224)
	i226     = fp_mul(i221,i225)
	i228     = fp_mul(fp_sqr(i226),i224)
	i229     = fp_mul(i225,i228)
	i230     = fp_mul(i222,i229)
	i232     = fp_mul(fp_sqr(i230),i223)
	i233     = fp_mul(i228,i232)
	i234     = fp_mul(i218,i233)
	i236     = fp_mul(fp_sqr(i234),i232)
	i237     = fp_mul(i233,i236)
	i238     = fp_mul(i230,i237)
	i239     = fp_mul(i234,i238)
	i240     = fp_mul(i226,i239)
	i241     = fp_mul(i236,i240)
	i244     = fp_mul(fp_mul(fp_sqr(i241),i241),i240)
	i245     = fp_mul(i239,i244)
	i246     = fp_mul(i244,i245)
	i248     = fp_mul(fp_sqr(i246),i237)
	i249     = fp_mul(i246,i248)
	i250     = fp_mul(i229,i249)
	i251     = fp_mul(i238,i250)
	i252     = fp_mul(i248,i251)
	i253     = fp_mul(i250,i252)
	i255     = fp_mul(fp_sqr(i253),i251)
	i256     = fp_mul(i252,i255)
	i257     = fp_mul(i245,i256)
	i258     = fp_mul(i253,i257)
	i260     = fp_mul(fp_sqr(i258),i255)
	i262     = fp_mul(fp_sqr(i260),i258)
	i263     = fp_mul(i256,i262)
	i264     = fp_mul(i260,i263)
	i265     = fp_mul(i257,i264)
	i266     = fp_mul(i262,i265)
	i267     = fp_mul(i249,i266)
	i268     = fp_mul(i263,i267)
	i269     = fp_mul(i266,i268)
	i270     = fp_mul(i241,i269)
	i271     = fp_mul(i269,i270)
	i272     = fp_mul(i267,i271)
	i273     = fp_mul(i265,i272)
	x123     = fp_mul(i264,i273)
	i275     = fp_mul(i271,x123)
	i276     = fp_mul(i270,i275)
	i277     = fp_mul(i268,i276)
	i280     = fp_mul(fp_mul(fp_sqr(i277),i277),i275)
	i281     = fp_mul(i272,i280)
	i282     = fp_mul(i273,i281)
	i283     = fp_mul(i281,i282)
	i284     = fp_mul(i280,i283)
	i285     = fp_mul(i276,i284)
	x128     = fp_mul(i277,i285)
	i676     = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i283,5444517870735015415413993718908291383296),i284),1),170141183460469231731687303715884105728),i282),1),680564733841876926926749214863536422912)
	i935     = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(i285,i676),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1321    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i935,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i1580    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1321),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i1966    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i1580,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2225    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i1966),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i2611    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2225,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i2870    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i2611),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i3256    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i2870,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i3515    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i3256),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	i3901    = fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(i3515,340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456)
	i4160    = fp_mul(fp_exp(fp_exp(fp_mul(fp_exp(fp_exp(fp_mul(x128,i3901),1),340282366920938463463374607431768211456),x128),1),340282366920938463463374607431768211456),x128)
	return      fp_mul(fp_exp(i4160,10633823966279326983230456482242756608),x123)
