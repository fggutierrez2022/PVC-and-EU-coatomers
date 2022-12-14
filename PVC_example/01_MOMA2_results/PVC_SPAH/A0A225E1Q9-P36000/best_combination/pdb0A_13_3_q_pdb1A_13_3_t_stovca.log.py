from pymol.cmd import *
code1 = 'pdb0A_13_3_q'
code2 = 'pdb1A_13_3_t'
cmd.load(code1 + '.pdb' , code1)
cmd.load(code2 + '.pdb' , code2)
cmd.hide('lines', 'all')
cmd.color('blue', code1)
cmd.color('green', code2)
cmd.show('sticks', 'all')
cmd.select('pair_1', '(' + code1 + ' and resi 573) + (' + code2 +' and resi 251)')
cmd.color('red', code1 + ' and resi 573')
cmd.color('orange',  code2 + ' and resi 251')
cmd.select('pair_2', '(' + code1 + ' and resi 574) + (' + code2 +' and resi 252)')
cmd.color('red', code1 + ' and resi 574')
cmd.color('orange',  code2 + ' and resi 252')
cmd.select('pair_3', '(' + code1 + ' and resi 575) + (' + code2 +' and resi 253)')
cmd.color('red', code1 + ' and resi 575')
cmd.color('orange',  code2 + ' and resi 253')
cmd.select('pair_4', '(' + code1 + ' and resi 576) + (' + code2 +' and resi 254)')
cmd.color('red', code1 + ' and resi 576')
cmd.color('orange',  code2 + ' and resi 254')
cmd.select('pair_5', '(' + code1 + ' and resi 577) + (' + code2 +' and resi 255)')
cmd.color('red', code1 + ' and resi 577')
cmd.color('orange',  code2 + ' and resi 255')
cmd.select('pair_6', '(' + code1 + ' and resi 578) + (' + code2 +' and resi 256)')
cmd.color('red', code1 + ' and resi 578')
cmd.color('orange',  code2 + ' and resi 256')
cmd.select('pair_7', '(' + code1 + ' and resi 579) + (' + code2 +' and resi 257)')
cmd.color('red', code1 + ' and resi 579')
cmd.color('orange',  code2 + ' and resi 257')
cmd.select('pair_8', '(' + code1 + ' and resi 580) + (' + code2 +' and resi 258)')
cmd.color('red', code1 + ' and resi 580')
cmd.color('orange',  code2 + ' and resi 258')
cmd.select('pair_9', '(' + code1 + ' and resi 581) + (' + code2 +' and resi 259)')
cmd.color('red', code1 + ' and resi 581')
cmd.color('orange',  code2 + ' and resi 259')
cmd.select('pair_10', '(' + code1 + ' and resi 582) + (' + code2 +' and resi 260)')
cmd.color('red', code1 + ' and resi 582')
cmd.color('orange',  code2 + ' and resi 260')
cmd.select('pair_11', '(' + code1 + ' and resi 583) + (' + code2 +' and resi 261)')
cmd.color('red', code1 + ' and resi 583')
cmd.color('orange',  code2 + ' and resi 261')
cmd.select('pair_12', '(' + code1 + ' and resi 584) + (' + code2 +' and resi 262)')
cmd.color('red', code1 + ' and resi 584')
cmd.color('orange',  code2 + ' and resi 262')
cmd.select('pair_13', '(' + code1 + ' and resi 585) + (' + code2 +' and resi 263)')
cmd.color('red', code1 + ' and resi 585')
cmd.color('orange',  code2 + ' and resi 263')
cmd.select('pair_14', '(' + code1 + ' and resi 586) + (' + code2 +' and resi 264)')
cmd.color('red', code1 + ' and resi 586')
cmd.color('orange',  code2 + ' and resi 264')
cmd.select('pair_15', '(' + code1 + ' and resi 587) + (' + code2 +' and resi 265)')
cmd.color('red', code1 + ' and resi 587')
cmd.color('orange',  code2 + ' and resi 265')
cmd.select('pair_16', '(' + code1 + ' and resi 590) + (' + code2 +' and resi 266)')
cmd.color('red', code1 + ' and resi 590')
cmd.color('orange',  code2 + ' and resi 266')
cmd.select('pair_17', '(' + code1 + ' and resi 591) + (' + code2 +' and resi 267)')
cmd.color('red', code1 + ' and resi 591')
cmd.color('orange',  code2 + ' and resi 267')
cmd.select('pair_18', '(' + code1 + ' and resi 592) + (' + code2 +' and resi 268)')
cmd.color('red', code1 + ' and resi 592')
cmd.color('orange',  code2 + ' and resi 268')
cmd.select('pair_19', '(' + code1 + ' and resi 593) + (' + code2 +' and resi 269)')
cmd.color('red', code1 + ' and resi 593')
cmd.color('orange',  code2 + ' and resi 269')
cmd.select('pair_20', '(' + code1 + ' and resi 594) + (' + code2 +' and resi 270)')
cmd.color('red', code1 + ' and resi 594')
cmd.color('orange',  code2 + ' and resi 270')
cmd.select('pair_21', '(' + code1 + ' and resi 595) + (' + code2 +' and resi 271)')
cmd.color('red', code1 + ' and resi 595')
cmd.color('orange',  code2 + ' and resi 271')
cmd.select('pair_22', '(' + code1 + ' and resi 596) + (' + code2 +' and resi 272)')
cmd.color('red', code1 + ' and resi 596')
cmd.color('orange',  code2 + ' and resi 272')
cmd.select('pair_23', '(' + code1 + ' and resi 597) + (' + code2 +' and resi 273)')
cmd.color('red', code1 + ' and resi 597')
cmd.color('orange',  code2 + ' and resi 273')
cmd.select('pair_24', '(' + code1 + ' and resi 598) + (' + code2 +' and resi 274)')
cmd.color('red', code1 + ' and resi 598')
cmd.color('orange',  code2 + ' and resi 274')
cmd.select('pair_25', '(' + code1 + ' and resi 599) + (' + code2 +' and resi 275)')
cmd.color('red', code1 + ' and resi 599')
cmd.color('orange',  code2 + ' and resi 275')
cmd.select('pair_26', '(' + code1 + ' and resi 600) + (' + code2 +' and resi 276)')
cmd.color('red', code1 + ' and resi 600')
cmd.color('orange',  code2 + ' and resi 276')
cmd.select('pair_27', '(' + code1 + ' and resi 601) + (' + code2 +' and resi 277)')
cmd.color('red', code1 + ' and resi 601')
cmd.color('orange',  code2 + ' and resi 277')
cmd.select('pair_28', '(' + code1 + ' and resi 602) + (' + code2 +' and resi 278)')
cmd.color('red', code1 + ' and resi 602')
cmd.color('orange',  code2 + ' and resi 278')
cmd.select('pair_29', '(' + code1 + ' and resi 603) + (' + code2 +' and resi 279)')
cmd.color('red', code1 + ' and resi 603')
cmd.color('orange',  code2 + ' and resi 279')
cmd.select('pair_30', '(' + code1 + ' and resi 604) + (' + code2 +' and resi 280)')
cmd.color('red', code1 + ' and resi 604')
cmd.color('orange',  code2 + ' and resi 280')
cmd.select('pair_31', '(' + code1 + ' and resi 605) + (' + code2 +' and resi 281)')
cmd.color('red', code1 + ' and resi 605')
cmd.color('orange',  code2 + ' and resi 281')
cmd.select('pair_32', '(' + code1 + ' and resi 606) + (' + code2 +' and resi 282)')
cmd.color('red', code1 + ' and resi 606')
cmd.color('orange',  code2 + ' and resi 282')
cmd.select('pair_33', '(' + code1 + ' and resi 607) + (' + code2 +' and resi 283)')
cmd.color('red', code1 + ' and resi 607')
cmd.color('orange',  code2 + ' and resi 283')
cmd.select('pair_34', '(' + code1 + ' and resi 638) + (' + code2 +' and resi 292)')
cmd.color('red', code1 + ' and resi 638')
cmd.color('orange',  code2 + ' and resi 292')
cmd.select('pair_35', '(' + code1 + ' and resi 639) + (' + code2 +' and resi 293)')
cmd.color('red', code1 + ' and resi 639')
cmd.color('orange',  code2 + ' and resi 293')
cmd.select('pair_36', '(' + code1 + ' and resi 640) + (' + code2 +' and resi 294)')
cmd.color('red', code1 + ' and resi 640')
cmd.color('orange',  code2 + ' and resi 294')
cmd.select('pair_37', '(' + code1 + ' and resi 641) + (' + code2 +' and resi 295)')
cmd.color('red', code1 + ' and resi 641')
cmd.color('orange',  code2 + ' and resi 295')
cmd.select('pair_38', '(' + code1 + ' and resi 642) + (' + code2 +' and resi 296)')
cmd.color('red', code1 + ' and resi 642')
cmd.color('orange',  code2 + ' and resi 296')
cmd.select('pair_39', '(' + code1 + ' and resi 643) + (' + code2 +' and resi 297)')
cmd.color('red', code1 + ' and resi 643')
cmd.color('orange',  code2 + ' and resi 297')
cmd.select('pair_40', '(' + code1 + ' and resi 644) + (' + code2 +' and resi 298)')
cmd.color('red', code1 + ' and resi 644')
cmd.color('orange',  code2 + ' and resi 298')
cmd.select('pair_41', '(' + code1 + ' and resi 645) + (' + code2 +' and resi 299)')
cmd.color('red', code1 + ' and resi 645')
cmd.color('orange',  code2 + ' and resi 299')
cmd.select('pair_42', '(' + code1 + ' and resi 646) + (' + code2 +' and resi 300)')
cmd.color('red', code1 + ' and resi 646')
cmd.color('orange',  code2 + ' and resi 300')
cmd.select('pair_43', '(' + code1 + ' and resi 647) + (' + code2 +' and resi 301)')
cmd.color('red', code1 + ' and resi 647')
cmd.color('orange',  code2 + ' and resi 301')
cmd.select('pair_44', '(' + code1 + ' and resi 648) + (' + code2 +' and resi 302)')
cmd.color('red', code1 + ' and resi 648')
cmd.color('orange',  code2 + ' and resi 302')
cmd.bg_color('gray70')
