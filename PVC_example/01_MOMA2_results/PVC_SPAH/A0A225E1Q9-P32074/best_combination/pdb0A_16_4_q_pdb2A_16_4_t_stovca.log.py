from pymol.cmd import *
code1 = 'pdb0A_16_4_q'
code2 = 'pdb2A_16_4_t'
cmd.load(code1 + '.pdb' , code1)
cmd.load(code2 + '.pdb' , code2)
cmd.hide('lines', 'all')
cmd.color('blue', code1)
cmd.color('green', code2)
cmd.show('sticks', 'all')
cmd.select('pair_1', '(' + code1 + ' and resi 650) + (' + code2 +' and resi 318)')
cmd.color('red', code1 + ' and resi 650')
cmd.color('orange',  code2 + ' and resi 318')
cmd.select('pair_2', '(' + code1 + ' and resi 651) + (' + code2 +' and resi 319)')
cmd.color('red', code1 + ' and resi 651')
cmd.color('orange',  code2 + ' and resi 319')
cmd.select('pair_3', '(' + code1 + ' and resi 652) + (' + code2 +' and resi 320)')
cmd.color('red', code1 + ' and resi 652')
cmd.color('orange',  code2 + ' and resi 320')
cmd.select('pair_4', '(' + code1 + ' and resi 653) + (' + code2 +' and resi 321)')
cmd.color('red', code1 + ' and resi 653')
cmd.color('orange',  code2 + ' and resi 321')
cmd.select('pair_5', '(' + code1 + ' and resi 654) + (' + code2 +' and resi 322)')
cmd.color('red', code1 + ' and resi 654')
cmd.color('orange',  code2 + ' and resi 322')
cmd.select('pair_6', '(' + code1 + ' and resi 655) + (' + code2 +' and resi 323)')
cmd.color('red', code1 + ' and resi 655')
cmd.color('orange',  code2 + ' and resi 323')
cmd.select('pair_7', '(' + code1 + ' and resi 656) + (' + code2 +' and resi 324)')
cmd.color('red', code1 + ' and resi 656')
cmd.color('orange',  code2 + ' and resi 324')
cmd.select('pair_8', '(' + code1 + ' and resi 657) + (' + code2 +' and resi 325)')
cmd.color('red', code1 + ' and resi 657')
cmd.color('orange',  code2 + ' and resi 325')
cmd.select('pair_9', '(' + code1 + ' and resi 658) + (' + code2 +' and resi 326)')
cmd.color('red', code1 + ' and resi 658')
cmd.color('orange',  code2 + ' and resi 326')
cmd.select('pair_10', '(' + code1 + ' and resi 659) + (' + code2 +' and resi 327)')
cmd.color('red', code1 + ' and resi 659')
cmd.color('orange',  code2 + ' and resi 327')
cmd.select('pair_11', '(' + code1 + ' and resi 660) + (' + code2 +' and resi 328)')
cmd.color('red', code1 + ' and resi 660')
cmd.color('orange',  code2 + ' and resi 328')
cmd.select('pair_12', '(' + code1 + ' and resi 661) + (' + code2 +' and resi 329)')
cmd.color('red', code1 + ' and resi 661')
cmd.color('orange',  code2 + ' and resi 329')
cmd.select('pair_13', '(' + code1 + ' and resi 662) + (' + code2 +' and resi 330)')
cmd.color('red', code1 + ' and resi 662')
cmd.color('orange',  code2 + ' and resi 330')
cmd.select('pair_14', '(' + code1 + ' and resi 663) + (' + code2 +' and resi 331)')
cmd.color('red', code1 + ' and resi 663')
cmd.color('orange',  code2 + ' and resi 331')
cmd.select('pair_15', '(' + code1 + ' and resi 664) + (' + code2 +' and resi 332)')
cmd.color('red', code1 + ' and resi 664')
cmd.color('orange',  code2 + ' and resi 332')
cmd.select('pair_16', '(' + code1 + ' and resi 665) + (' + code2 +' and resi 333)')
cmd.color('red', code1 + ' and resi 665')
cmd.color('orange',  code2 + ' and resi 333')
cmd.select('pair_17', '(' + code1 + ' and resi 666) + (' + code2 +' and resi 340)')
cmd.color('red', code1 + ' and resi 666')
cmd.color('orange',  code2 + ' and resi 340')
cmd.select('pair_18', '(' + code1 + ' and resi 667) + (' + code2 +' and resi 341)')
cmd.color('red', code1 + ' and resi 667')
cmd.color('orange',  code2 + ' and resi 341')
cmd.select('pair_19', '(' + code1 + ' and resi 668) + (' + code2 +' and resi 342)')
cmd.color('red', code1 + ' and resi 668')
cmd.color('orange',  code2 + ' and resi 342')
cmd.select('pair_20', '(' + code1 + ' and resi 669) + (' + code2 +' and resi 343)')
cmd.color('red', code1 + ' and resi 669')
cmd.color('orange',  code2 + ' and resi 343')
cmd.select('pair_21', '(' + code1 + ' and resi 670) + (' + code2 +' and resi 344)')
cmd.color('red', code1 + ' and resi 670')
cmd.color('orange',  code2 + ' and resi 344')
cmd.select('pair_22', '(' + code1 + ' and resi 672) + (' + code2 +' and resi 346)')
cmd.color('red', code1 + ' and resi 672')
cmd.color('orange',  code2 + ' and resi 346')
cmd.select('pair_23', '(' + code1 + ' and resi 673) + (' + code2 +' and resi 347)')
cmd.color('red', code1 + ' and resi 673')
cmd.color('orange',  code2 + ' and resi 347')
cmd.select('pair_24', '(' + code1 + ' and resi 674) + (' + code2 +' and resi 348)')
cmd.color('red', code1 + ' and resi 674')
cmd.color('orange',  code2 + ' and resi 348')
cmd.select('pair_25', '(' + code1 + ' and resi 675) + (' + code2 +' and resi 349)')
cmd.color('red', code1 + ' and resi 675')
cmd.color('orange',  code2 + ' and resi 349')
cmd.select('pair_26', '(' + code1 + ' and resi 676) + (' + code2 +' and resi 350)')
cmd.color('red', code1 + ' and resi 676')
cmd.color('orange',  code2 + ' and resi 350')
cmd.select('pair_27', '(' + code1 + ' and resi 677) + (' + code2 +' and resi 351)')
cmd.color('red', code1 + ' and resi 677')
cmd.color('orange',  code2 + ' and resi 351')
cmd.select('pair_28', '(' + code1 + ' and resi 678) + (' + code2 +' and resi 352)')
cmd.color('red', code1 + ' and resi 678')
cmd.color('orange',  code2 + ' and resi 352')
cmd.select('pair_29', '(' + code1 + ' and resi 680) + (' + code2 +' and resi 353)')
cmd.color('red', code1 + ' and resi 680')
cmd.color('orange',  code2 + ' and resi 353')
cmd.select('pair_30', '(' + code1 + ' and resi 681) + (' + code2 +' and resi 354)')
cmd.color('red', code1 + ' and resi 681')
cmd.color('orange',  code2 + ' and resi 354')
cmd.select('pair_31', '(' + code1 + ' and resi 682) + (' + code2 +' and resi 355)')
cmd.color('red', code1 + ' and resi 682')
cmd.color('orange',  code2 + ' and resi 355')
cmd.select('pair_32', '(' + code1 + ' and resi 683) + (' + code2 +' and resi 356)')
cmd.color('red', code1 + ' and resi 683')
cmd.color('orange',  code2 + ' and resi 356')
cmd.select('pair_33', '(' + code1 + ' and resi 684) + (' + code2 +' and resi 357)')
cmd.color('red', code1 + ' and resi 684')
cmd.color('orange',  code2 + ' and resi 357')
cmd.select('pair_34', '(' + code1 + ' and resi 685) + (' + code2 +' and resi 358)')
cmd.color('red', code1 + ' and resi 685')
cmd.color('orange',  code2 + ' and resi 358')
cmd.select('pair_35', '(' + code1 + ' and resi 686) + (' + code2 +' and resi 359)')
cmd.color('red', code1 + ' and resi 686')
cmd.color('orange',  code2 + ' and resi 359')
cmd.select('pair_36', '(' + code1 + ' and resi 687) + (' + code2 +' and resi 360)')
cmd.color('red', code1 + ' and resi 687')
cmd.color('orange',  code2 + ' and resi 360')
cmd.select('pair_37', '(' + code1 + ' and resi 688) + (' + code2 +' and resi 361)')
cmd.color('red', code1 + ' and resi 688')
cmd.color('orange',  code2 + ' and resi 361')
cmd.select('pair_38', '(' + code1 + ' and resi 689) + (' + code2 +' and resi 362)')
cmd.color('red', code1 + ' and resi 689')
cmd.color('orange',  code2 + ' and resi 362')
cmd.select('pair_39', '(' + code1 + ' and resi 690) + (' + code2 +' and resi 363)')
cmd.color('red', code1 + ' and resi 690')
cmd.color('orange',  code2 + ' and resi 363')
cmd.select('pair_40', '(' + code1 + ' and resi 691) + (' + code2 +' and resi 364)')
cmd.color('red', code1 + ' and resi 691')
cmd.color('orange',  code2 + ' and resi 364')
cmd.select('pair_41', '(' + code1 + ' and resi 692) + (' + code2 +' and resi 365)')
cmd.color('red', code1 + ' and resi 692')
cmd.color('orange',  code2 + ' and resi 365')
cmd.select('pair_42', '(' + code1 + ' and resi 693) + (' + code2 +' and resi 366)')
cmd.color('red', code1 + ' and resi 693')
cmd.color('orange',  code2 + ' and resi 366')
cmd.select('pair_43', '(' + code1 + ' and resi 694) + (' + code2 +' and resi 367)')
cmd.color('red', code1 + ' and resi 694')
cmd.color('orange',  code2 + ' and resi 367')
cmd.select('pair_44', '(' + code1 + ' and resi 695) + (' + code2 +' and resi 368)')
cmd.color('red', code1 + ' and resi 695')
cmd.color('orange',  code2 + ' and resi 368')
cmd.select('pair_45', '(' + code1 + ' and resi 696) + (' + code2 +' and resi 369)')
cmd.color('red', code1 + ' and resi 696')
cmd.color('orange',  code2 + ' and resi 369')
cmd.select('pair_46', '(' + code1 + ' and resi 702) + (' + code2 +' and resi 376)')
cmd.color('red', code1 + ' and resi 702')
cmd.color('orange',  code2 + ' and resi 376')
cmd.select('pair_47', '(' + code1 + ' and resi 703) + (' + code2 +' and resi 377)')
cmd.color('red', code1 + ' and resi 703')
cmd.color('orange',  code2 + ' and resi 377')
cmd.select('pair_48', '(' + code1 + ' and resi 704) + (' + code2 +' and resi 378)')
cmd.color('red', code1 + ' and resi 704')
cmd.color('orange',  code2 + ' and resi 378')
cmd.select('pair_49', '(' + code1 + ' and resi 705) + (' + code2 +' and resi 379)')
cmd.color('red', code1 + ' and resi 705')
cmd.color('orange',  code2 + ' and resi 379')
cmd.select('pair_50', '(' + code1 + ' and resi 706) + (' + code2 +' and resi 380)')
cmd.color('red', code1 + ' and resi 706')
cmd.color('orange',  code2 + ' and resi 380')
cmd.select('pair_51', '(' + code1 + ' and resi 707) + (' + code2 +' and resi 381)')
cmd.color('red', code1 + ' and resi 707')
cmd.color('orange',  code2 + ' and resi 381')
cmd.select('pair_52', '(' + code1 + ' and resi 708) + (' + code2 +' and resi 382)')
cmd.color('red', code1 + ' and resi 708')
cmd.color('orange',  code2 + ' and resi 382')
cmd.select('pair_53', '(' + code1 + ' and resi 709) + (' + code2 +' and resi 383)')
cmd.color('red', code1 + ' and resi 709')
cmd.color('orange',  code2 + ' and resi 383')
cmd.bg_color('gray70')
