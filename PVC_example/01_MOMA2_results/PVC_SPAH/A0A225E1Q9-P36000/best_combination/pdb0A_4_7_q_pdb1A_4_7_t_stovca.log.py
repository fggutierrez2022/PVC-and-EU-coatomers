from pymol.cmd import *
code1 = 'pdb0A_4_7_q'
code2 = 'pdb1A_4_7_t'
cmd.load(code1 + '.pdb' , code1)
cmd.load(code2 + '.pdb' , code2)
cmd.hide('lines', 'all')
cmd.color('blue', code1)
cmd.color('green', code2)
cmd.show('sticks', 'all')
cmd.select('pair_1', '(' + code1 + ' and resi 430) + (' + code2 +' and resi 75)')
cmd.color('red', code1 + ' and resi 430')
cmd.color('orange',  code2 + ' and resi 75')
cmd.select('pair_2', '(' + code1 + ' and resi 431) + (' + code2 +' and resi 76)')
cmd.color('red', code1 + ' and resi 431')
cmd.color('orange',  code2 + ' and resi 76')
cmd.select('pair_3', '(' + code1 + ' and resi 432) + (' + code2 +' and resi 77)')
cmd.color('red', code1 + ' and resi 432')
cmd.color('orange',  code2 + ' and resi 77')
cmd.select('pair_4', '(' + code1 + ' and resi 433) + (' + code2 +' and resi 78)')
cmd.color('red', code1 + ' and resi 433')
cmd.color('orange',  code2 + ' and resi 78')
cmd.select('pair_5', '(' + code1 + ' and resi 434) + (' + code2 +' and resi 79)')
cmd.color('red', code1 + ' and resi 434')
cmd.color('orange',  code2 + ' and resi 79')
cmd.select('pair_6', '(' + code1 + ' and resi 435) + (' + code2 +' and resi 80)')
cmd.color('red', code1 + ' and resi 435')
cmd.color('orange',  code2 + ' and resi 80')
cmd.select('pair_7', '(' + code1 + ' and resi 436) + (' + code2 +' and resi 81)')
cmd.color('red', code1 + ' and resi 436')
cmd.color('orange',  code2 + ' and resi 81')
cmd.select('pair_8', '(' + code1 + ' and resi 437) + (' + code2 +' and resi 82)')
cmd.color('red', code1 + ' and resi 437')
cmd.color('orange',  code2 + ' and resi 82')
cmd.select('pair_9', '(' + code1 + ' and resi 438) + (' + code2 +' and resi 83)')
cmd.color('red', code1 + ' and resi 438')
cmd.color('orange',  code2 + ' and resi 83')
cmd.select('pair_10', '(' + code1 + ' and resi 439) + (' + code2 +' and resi 84)')
cmd.color('red', code1 + ' and resi 439')
cmd.color('orange',  code2 + ' and resi 84')
cmd.select('pair_11', '(' + code1 + ' and resi 440) + (' + code2 +' and resi 85)')
cmd.color('red', code1 + ' and resi 440')
cmd.color('orange',  code2 + ' and resi 85')
cmd.select('pair_12', '(' + code1 + ' and resi 441) + (' + code2 +' and resi 86)')
cmd.color('red', code1 + ' and resi 441')
cmd.color('orange',  code2 + ' and resi 86')
cmd.select('pair_13', '(' + code1 + ' and resi 447) + (' + code2 +' and resi 97)')
cmd.color('red', code1 + ' and resi 447')
cmd.color('orange',  code2 + ' and resi 97')
cmd.select('pair_14', '(' + code1 + ' and resi 448) + (' + code2 +' and resi 98)')
cmd.color('red', code1 + ' and resi 448')
cmd.color('orange',  code2 + ' and resi 98')
cmd.select('pair_15', '(' + code1 + ' and resi 449) + (' + code2 +' and resi 99)')
cmd.color('red', code1 + ' and resi 449')
cmd.color('orange',  code2 + ' and resi 99')
cmd.select('pair_16', '(' + code1 + ' and resi 450) + (' + code2 +' and resi 100)')
cmd.color('red', code1 + ' and resi 450')
cmd.color('orange',  code2 + ' and resi 100')
cmd.select('pair_17', '(' + code1 + ' and resi 451) + (' + code2 +' and resi 101)')
cmd.color('red', code1 + ' and resi 451')
cmd.color('orange',  code2 + ' and resi 101')
cmd.select('pair_18', '(' + code1 + ' and resi 452) + (' + code2 +' and resi 102)')
cmd.color('red', code1 + ' and resi 452')
cmd.color('orange',  code2 + ' and resi 102')
cmd.select('pair_19', '(' + code1 + ' and resi 453) + (' + code2 +' and resi 103)')
cmd.color('red', code1 + ' and resi 453')
cmd.color('orange',  code2 + ' and resi 103')
cmd.select('pair_20', '(' + code1 + ' and resi 454) + (' + code2 +' and resi 104)')
cmd.color('red', code1 + ' and resi 454')
cmd.color('orange',  code2 + ' and resi 104')
cmd.select('pair_21', '(' + code1 + ' and resi 455) + (' + code2 +' and resi 105)')
cmd.color('red', code1 + ' and resi 455')
cmd.color('orange',  code2 + ' and resi 105')
cmd.select('pair_22', '(' + code1 + ' and resi 456) + (' + code2 +' and resi 106)')
cmd.color('red', code1 + ' and resi 456')
cmd.color('orange',  code2 + ' and resi 106')
cmd.select('pair_23', '(' + code1 + ' and resi 457) + (' + code2 +' and resi 107)')
cmd.color('red', code1 + ' and resi 457')
cmd.color('orange',  code2 + ' and resi 107')
cmd.select('pair_24', '(' + code1 + ' and resi 458) + (' + code2 +' and resi 108)')
cmd.color('red', code1 + ' and resi 458')
cmd.color('orange',  code2 + ' and resi 108')
cmd.select('pair_25', '(' + code1 + ' and resi 459) + (' + code2 +' and resi 109)')
cmd.color('red', code1 + ' and resi 459')
cmd.color('orange',  code2 + ' and resi 109')
cmd.select('pair_26', '(' + code1 + ' and resi 460) + (' + code2 +' and resi 110)')
cmd.color('red', code1 + ' and resi 460')
cmd.color('orange',  code2 + ' and resi 110')
cmd.select('pair_27', '(' + code1 + ' and resi 461) + (' + code2 +' and resi 111)')
cmd.color('red', code1 + ' and resi 461')
cmd.color('orange',  code2 + ' and resi 111')
cmd.select('pair_28', '(' + code1 + ' and resi 462) + (' + code2 +' and resi 112)')
cmd.color('red', code1 + ' and resi 462')
cmd.color('orange',  code2 + ' and resi 112')
cmd.select('pair_29', '(' + code1 + ' and resi 463) + (' + code2 +' and resi 113)')
cmd.color('red', code1 + ' and resi 463')
cmd.color('orange',  code2 + ' and resi 113')
cmd.select('pair_30', '(' + code1 + ' and resi 464) + (' + code2 +' and resi 114)')
cmd.color('red', code1 + ' and resi 464')
cmd.color('orange',  code2 + ' and resi 114')
cmd.select('pair_31', '(' + code1 + ' and resi 465) + (' + code2 +' and resi 115)')
cmd.color('red', code1 + ' and resi 465')
cmd.color('orange',  code2 + ' and resi 115')
cmd.select('pair_32', '(' + code1 + ' and resi 466) + (' + code2 +' and resi 116)')
cmd.color('red', code1 + ' and resi 466')
cmd.color('orange',  code2 + ' and resi 116')
cmd.select('pair_33', '(' + code1 + ' and resi 467) + (' + code2 +' and resi 117)')
cmd.color('red', code1 + ' and resi 467')
cmd.color('orange',  code2 + ' and resi 117')
cmd.select('pair_34', '(' + code1 + ' and resi 468) + (' + code2 +' and resi 118)')
cmd.color('red', code1 + ' and resi 468')
cmd.color('orange',  code2 + ' and resi 118')
cmd.select('pair_35', '(' + code1 + ' and resi 469) + (' + code2 +' and resi 119)')
cmd.color('red', code1 + ' and resi 469')
cmd.color('orange',  code2 + ' and resi 119')
cmd.select('pair_36', '(' + code1 + ' and resi 470) + (' + code2 +' and resi 120)')
cmd.color('red', code1 + ' and resi 470')
cmd.color('orange',  code2 + ' and resi 120')
cmd.select('pair_37', '(' + code1 + ' and resi 471) + (' + code2 +' and resi 121)')
cmd.color('red', code1 + ' and resi 471')
cmd.color('orange',  code2 + ' and resi 121')
cmd.select('pair_38', '(' + code1 + ' and resi 472) + (' + code2 +' and resi 122)')
cmd.color('red', code1 + ' and resi 472')
cmd.color('orange',  code2 + ' and resi 122')
cmd.select('pair_39', '(' + code1 + ' and resi 473) + (' + code2 +' and resi 123)')
cmd.color('red', code1 + ' and resi 473')
cmd.color('orange',  code2 + ' and resi 123')
cmd.select('pair_40', '(' + code1 + ' and resi 474) + (' + code2 +' and resi 124)')
cmd.color('red', code1 + ' and resi 474')
cmd.color('orange',  code2 + ' and resi 124')
cmd.select('pair_41', '(' + code1 + ' and resi 475) + (' + code2 +' and resi 125)')
cmd.color('red', code1 + ' and resi 475')
cmd.color('orange',  code2 + ' and resi 125')
cmd.select('pair_42', '(' + code1 + ' and resi 480) + (' + code2 +' and resi 134)')
cmd.color('red', code1 + ' and resi 480')
cmd.color('orange',  code2 + ' and resi 134')
cmd.select('pair_43', '(' + code1 + ' and resi 481) + (' + code2 +' and resi 135)')
cmd.color('red', code1 + ' and resi 481')
cmd.color('orange',  code2 + ' and resi 135')
cmd.select('pair_44', '(' + code1 + ' and resi 482) + (' + code2 +' and resi 136)')
cmd.color('red', code1 + ' and resi 482')
cmd.color('orange',  code2 + ' and resi 136')
cmd.select('pair_45', '(' + code1 + ' and resi 483) + (' + code2 +' and resi 137)')
cmd.color('red', code1 + ' and resi 483')
cmd.color('orange',  code2 + ' and resi 137')
cmd.select('pair_46', '(' + code1 + ' and resi 484) + (' + code2 +' and resi 138)')
cmd.color('red', code1 + ' and resi 484')
cmd.color('orange',  code2 + ' and resi 138')
cmd.select('pair_47', '(' + code1 + ' and resi 485) + (' + code2 +' and resi 139)')
cmd.color('red', code1 + ' and resi 485')
cmd.color('orange',  code2 + ' and resi 139')
cmd.select('pair_48', '(' + code1 + ' and resi 486) + (' + code2 +' and resi 140)')
cmd.color('red', code1 + ' and resi 486')
cmd.color('orange',  code2 + ' and resi 140')
cmd.select('pair_49', '(' + code1 + ' and resi 487) + (' + code2 +' and resi 141)')
cmd.color('red', code1 + ' and resi 487')
cmd.color('orange',  code2 + ' and resi 141')
cmd.select('pair_50', '(' + code1 + ' and resi 488) + (' + code2 +' and resi 142)')
cmd.color('red', code1 + ' and resi 488')
cmd.color('orange',  code2 + ' and resi 142')
cmd.select('pair_51', '(' + code1 + ' and resi 489) + (' + code2 +' and resi 143)')
cmd.color('red', code1 + ' and resi 489')
cmd.color('orange',  code2 + ' and resi 143')
cmd.select('pair_52', '(' + code1 + ' and resi 490) + (' + code2 +' and resi 144)')
cmd.color('red', code1 + ' and resi 490')
cmd.color('orange',  code2 + ' and resi 144')
cmd.select('pair_53', '(' + code1 + ' and resi 491) + (' + code2 +' and resi 145)')
cmd.color('red', code1 + ' and resi 491')
cmd.color('orange',  code2 + ' and resi 145')
cmd.select('pair_54', '(' + code1 + ' and resi 492) + (' + code2 +' and resi 146)')
cmd.color('red', code1 + ' and resi 492')
cmd.color('orange',  code2 + ' and resi 146')
cmd.select('pair_55', '(' + code1 + ' and resi 493) + (' + code2 +' and resi 147)')
cmd.color('red', code1 + ' and resi 493')
cmd.color('orange',  code2 + ' and resi 147')
cmd.select('pair_56', '(' + code1 + ' and resi 494) + (' + code2 +' and resi 148)')
cmd.color('red', code1 + ' and resi 494')
cmd.color('orange',  code2 + ' and resi 148')
cmd.select('pair_57', '(' + code1 + ' and resi 495) + (' + code2 +' and resi 149)')
cmd.color('red', code1 + ' and resi 495')
cmd.color('orange',  code2 + ' and resi 149')
cmd.select('pair_58', '(' + code1 + ' and resi 496) + (' + code2 +' and resi 150)')
cmd.color('red', code1 + ' and resi 496')
cmd.color('orange',  code2 + ' and resi 150')
cmd.select('pair_59', '(' + code1 + ' and resi 497) + (' + code2 +' and resi 151)')
cmd.color('red', code1 + ' and resi 497')
cmd.color('orange',  code2 + ' and resi 151')
cmd.select('pair_60', '(' + code1 + ' and resi 498) + (' + code2 +' and resi 152)')
cmd.color('red', code1 + ' and resi 498')
cmd.color('orange',  code2 + ' and resi 152')
cmd.select('pair_61', '(' + code1 + ' and resi 499) + (' + code2 +' and resi 153)')
cmd.color('red', code1 + ' and resi 499')
cmd.color('orange',  code2 + ' and resi 153')
cmd.select('pair_62', '(' + code1 + ' and resi 500) + (' + code2 +' and resi 154)')
cmd.color('red', code1 + ' and resi 500')
cmd.color('orange',  code2 + ' and resi 154')
cmd.select('pair_63', '(' + code1 + ' and resi 501) + (' + code2 +' and resi 155)')
cmd.color('red', code1 + ' and resi 501')
cmd.color('orange',  code2 + ' and resi 155')
cmd.select('pair_64', '(' + code1 + ' and resi 502) + (' + code2 +' and resi 156)')
cmd.color('red', code1 + ' and resi 502')
cmd.color('orange',  code2 + ' and resi 156')
cmd.select('pair_65', '(' + code1 + ' and resi 503) + (' + code2 +' and resi 157)')
cmd.color('red', code1 + ' and resi 503')
cmd.color('orange',  code2 + ' and resi 157')
cmd.select('pair_66', '(' + code1 + ' and resi 504) + (' + code2 +' and resi 158)')
cmd.color('red', code1 + ' and resi 504')
cmd.color('orange',  code2 + ' and resi 158')
cmd.select('pair_67', '(' + code1 + ' and resi 505) + (' + code2 +' and resi 159)')
cmd.color('red', code1 + ' and resi 505')
cmd.color('orange',  code2 + ' and resi 159')
cmd.select('pair_68', '(' + code1 + ' and resi 506) + (' + code2 +' and resi 160)')
cmd.color('red', code1 + ' and resi 506')
cmd.color('orange',  code2 + ' and resi 160')
cmd.select('pair_69', '(' + code1 + ' and resi 507) + (' + code2 +' and resi 161)')
cmd.color('red', code1 + ' and resi 507')
cmd.color('orange',  code2 + ' and resi 161')
cmd.select('pair_70', '(' + code1 + ' and resi 508) + (' + code2 +' and resi 162)')
cmd.color('red', code1 + ' and resi 508')
cmd.color('orange',  code2 + ' and resi 162')
cmd.select('pair_71', '(' + code1 + ' and resi 511) + (' + code2 +' and resi 173)')
cmd.color('red', code1 + ' and resi 511')
cmd.color('orange',  code2 + ' and resi 173')
cmd.select('pair_72', '(' + code1 + ' and resi 512) + (' + code2 +' and resi 174)')
cmd.color('red', code1 + ' and resi 512')
cmd.color('orange',  code2 + ' and resi 174')
cmd.select('pair_73', '(' + code1 + ' and resi 513) + (' + code2 +' and resi 175)')
cmd.color('red', code1 + ' and resi 513')
cmd.color('orange',  code2 + ' and resi 175')
cmd.select('pair_74', '(' + code1 + ' and resi 514) + (' + code2 +' and resi 176)')
cmd.color('red', code1 + ' and resi 514')
cmd.color('orange',  code2 + ' and resi 176')
cmd.select('pair_75', '(' + code1 + ' and resi 515) + (' + code2 +' and resi 177)')
cmd.color('red', code1 + ' and resi 515')
cmd.color('orange',  code2 + ' and resi 177')
cmd.select('pair_76', '(' + code1 + ' and resi 516) + (' + code2 +' and resi 178)')
cmd.color('red', code1 + ' and resi 516')
cmd.color('orange',  code2 + ' and resi 178')
cmd.select('pair_77', '(' + code1 + ' and resi 517) + (' + code2 +' and resi 179)')
cmd.color('red', code1 + ' and resi 517')
cmd.color('orange',  code2 + ' and resi 179')
cmd.select('pair_78', '(' + code1 + ' and resi 518) + (' + code2 +' and resi 180)')
cmd.color('red', code1 + ' and resi 518')
cmd.color('orange',  code2 + ' and resi 180')
cmd.select('pair_79', '(' + code1 + ' and resi 519) + (' + code2 +' and resi 181)')
cmd.color('red', code1 + ' and resi 519')
cmd.color('orange',  code2 + ' and resi 181')
cmd.select('pair_80', '(' + code1 + ' and resi 520) + (' + code2 +' and resi 182)')
cmd.color('red', code1 + ' and resi 520')
cmd.color('orange',  code2 + ' and resi 182')
cmd.select('pair_81', '(' + code1 + ' and resi 521) + (' + code2 +' and resi 183)')
cmd.color('red', code1 + ' and resi 521')
cmd.color('orange',  code2 + ' and resi 183')
cmd.select('pair_82', '(' + code1 + ' and resi 522) + (' + code2 +' and resi 184)')
cmd.color('red', code1 + ' and resi 522')
cmd.color('orange',  code2 + ' and resi 184')
cmd.select('pair_83', '(' + code1 + ' and resi 523) + (' + code2 +' and resi 185)')
cmd.color('red', code1 + ' and resi 523')
cmd.color('orange',  code2 + ' and resi 185')
cmd.select('pair_84', '(' + code1 + ' and resi 524) + (' + code2 +' and resi 186)')
cmd.color('red', code1 + ' and resi 524')
cmd.color('orange',  code2 + ' and resi 186')
cmd.select('pair_85', '(' + code1 + ' and resi 525) + (' + code2 +' and resi 187)')
cmd.color('red', code1 + ' and resi 525')
cmd.color('orange',  code2 + ' and resi 187')
cmd.select('pair_86', '(' + code1 + ' and resi 526) + (' + code2 +' and resi 188)')
cmd.color('red', code1 + ' and resi 526')
cmd.color('orange',  code2 + ' and resi 188')
cmd.select('pair_87', '(' + code1 + ' and resi 527) + (' + code2 +' and resi 189)')
cmd.color('red', code1 + ' and resi 527')
cmd.color('orange',  code2 + ' and resi 189')
cmd.select('pair_88', '(' + code1 + ' and resi 528) + (' + code2 +' and resi 190)')
cmd.color('red', code1 + ' and resi 528')
cmd.color('orange',  code2 + ' and resi 190')
cmd.select('pair_89', '(' + code1 + ' and resi 529) + (' + code2 +' and resi 191)')
cmd.color('red', code1 + ' and resi 529')
cmd.color('orange',  code2 + ' and resi 191')
cmd.select('pair_90', '(' + code1 + ' and resi 530) + (' + code2 +' and resi 192)')
cmd.color('red', code1 + ' and resi 530')
cmd.color('orange',  code2 + ' and resi 192')
cmd.select('pair_91', '(' + code1 + ' and resi 531) + (' + code2 +' and resi 193)')
cmd.color('red', code1 + ' and resi 531')
cmd.color('orange',  code2 + ' and resi 193')
cmd.select('pair_92', '(' + code1 + ' and resi 532) + (' + code2 +' and resi 194)')
cmd.color('red', code1 + ' and resi 532')
cmd.color('orange',  code2 + ' and resi 194')
cmd.select('pair_93', '(' + code1 + ' and resi 533) + (' + code2 +' and resi 195)')
cmd.color('red', code1 + ' and resi 533')
cmd.color('orange',  code2 + ' and resi 195')
cmd.select('pair_94', '(' + code1 + ' and resi 534) + (' + code2 +' and resi 196)')
cmd.color('red', code1 + ' and resi 534')
cmd.color('orange',  code2 + ' and resi 196')
cmd.select('pair_95', '(' + code1 + ' and resi 535) + (' + code2 +' and resi 197)')
cmd.color('red', code1 + ' and resi 535')
cmd.color('orange',  code2 + ' and resi 197')
cmd.select('pair_96', '(' + code1 + ' and resi 536) + (' + code2 +' and resi 198)')
cmd.color('red', code1 + ' and resi 536')
cmd.color('orange',  code2 + ' and resi 198')
cmd.select('pair_97', '(' + code1 + ' and resi 537) + (' + code2 +' and resi 199)')
cmd.color('red', code1 + ' and resi 537')
cmd.color('orange',  code2 + ' and resi 199')
cmd.select('pair_98', '(' + code1 + ' and resi 538) + (' + code2 +' and resi 200)')
cmd.color('red', code1 + ' and resi 538')
cmd.color('orange',  code2 + ' and resi 200')
cmd.bg_color('gray70')