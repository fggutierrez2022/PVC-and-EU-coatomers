from pymol.cmd import *
code1 = 'pdb0A_16_7_q'
code2 = 'pdb2A_16_7_t'
cmd.load(code1 + '.pdb' , code1)
cmd.load(code2 + '.pdb' , code2)
cmd.hide('lines', 'all')
cmd.color('blue', code1)
cmd.color('green', code2)
cmd.show('sticks', 'all')
cmd.select('pair_1', '(' + code1 + ' and resi 268) + (' + code2 +' and resi 306)')
cmd.color('red', code1 + ' and resi 268')
cmd.color('orange',  code2 + ' and resi 306')
cmd.select('pair_2', '(' + code1 + ' and resi 269) + (' + code2 +' and resi 307)')
cmd.color('red', code1 + ' and resi 269')
cmd.color('orange',  code2 + ' and resi 307')
cmd.select('pair_3', '(' + code1 + ' and resi 270) + (' + code2 +' and resi 308)')
cmd.color('red', code1 + ' and resi 270')
cmd.color('orange',  code2 + ' and resi 308')
cmd.select('pair_4', '(' + code1 + ' and resi 271) + (' + code2 +' and resi 309)')
cmd.color('red', code1 + ' and resi 271')
cmd.color('orange',  code2 + ' and resi 309')
cmd.select('pair_5', '(' + code1 + ' and resi 272) + (' + code2 +' and resi 310)')
cmd.color('red', code1 + ' and resi 272')
cmd.color('orange',  code2 + ' and resi 310')
cmd.select('pair_6', '(' + code1 + ' and resi 273) + (' + code2 +' and resi 311)')
cmd.color('red', code1 + ' and resi 273')
cmd.color('orange',  code2 + ' and resi 311')
cmd.select('pair_7', '(' + code1 + ' and resi 274) + (' + code2 +' and resi 312)')
cmd.color('red', code1 + ' and resi 274')
cmd.color('orange',  code2 + ' and resi 312')
cmd.select('pair_8', '(' + code1 + ' and resi 284) + (' + code2 +' and resi 316)')
cmd.color('red', code1 + ' and resi 284')
cmd.color('orange',  code2 + ' and resi 316')
cmd.select('pair_9', '(' + code1 + ' and resi 285) + (' + code2 +' and resi 317)')
cmd.color('red', code1 + ' and resi 285')
cmd.color('orange',  code2 + ' and resi 317')
cmd.select('pair_10', '(' + code1 + ' and resi 286) + (' + code2 +' and resi 318)')
cmd.color('red', code1 + ' and resi 286')
cmd.color('orange',  code2 + ' and resi 318')
cmd.select('pair_11', '(' + code1 + ' and resi 287) + (' + code2 +' and resi 319)')
cmd.color('red', code1 + ' and resi 287')
cmd.color('orange',  code2 + ' and resi 319')
cmd.select('pair_12', '(' + code1 + ' and resi 288) + (' + code2 +' and resi 320)')
cmd.color('red', code1 + ' and resi 288')
cmd.color('orange',  code2 + ' and resi 320')
cmd.select('pair_13', '(' + code1 + ' and resi 289) + (' + code2 +' and resi 321)')
cmd.color('red', code1 + ' and resi 289')
cmd.color('orange',  code2 + ' and resi 321')
cmd.select('pair_14', '(' + code1 + ' and resi 294) + (' + code2 +' and resi 325)')
cmd.color('red', code1 + ' and resi 294')
cmd.color('orange',  code2 + ' and resi 325')
cmd.select('pair_15', '(' + code1 + ' and resi 295) + (' + code2 +' and resi 326)')
cmd.color('red', code1 + ' and resi 295')
cmd.color('orange',  code2 + ' and resi 326')
cmd.select('pair_16', '(' + code1 + ' and resi 296) + (' + code2 +' and resi 327)')
cmd.color('red', code1 + ' and resi 296')
cmd.color('orange',  code2 + ' and resi 327')
cmd.select('pair_17', '(' + code1 + ' and resi 297) + (' + code2 +' and resi 328)')
cmd.color('red', code1 + ' and resi 297')
cmd.color('orange',  code2 + ' and resi 328')
cmd.select('pair_18', '(' + code1 + ' and resi 298) + (' + code2 +' and resi 329)')
cmd.color('red', code1 + ' and resi 298')
cmd.color('orange',  code2 + ' and resi 329')
cmd.select('pair_19', '(' + code1 + ' and resi 299) + (' + code2 +' and resi 330)')
cmd.color('red', code1 + ' and resi 299')
cmd.color('orange',  code2 + ' and resi 330')
cmd.select('pair_20', '(' + code1 + ' and resi 300) + (' + code2 +' and resi 331)')
cmd.color('red', code1 + ' and resi 300')
cmd.color('orange',  code2 + ' and resi 331')
cmd.select('pair_21', '(' + code1 + ' and resi 301) + (' + code2 +' and resi 332)')
cmd.color('red', code1 + ' and resi 301')
cmd.color('orange',  code2 + ' and resi 332')
cmd.select('pair_22', '(' + code1 + ' and resi 302) + (' + code2 +' and resi 333)')
cmd.color('red', code1 + ' and resi 302')
cmd.color('orange',  code2 + ' and resi 333')
cmd.select('pair_23', '(' + code1 + ' and resi 307) + (' + code2 +' and resi 354)')
cmd.color('red', code1 + ' and resi 307')
cmd.color('orange',  code2 + ' and resi 354')
cmd.select('pair_24', '(' + code1 + ' and resi 308) + (' + code2 +' and resi 355)')
cmd.color('red', code1 + ' and resi 308')
cmd.color('orange',  code2 + ' and resi 355')
cmd.select('pair_25', '(' + code1 + ' and resi 309) + (' + code2 +' and resi 356)')
cmd.color('red', code1 + ' and resi 309')
cmd.color('orange',  code2 + ' and resi 356')
cmd.select('pair_26', '(' + code1 + ' and resi 310) + (' + code2 +' and resi 357)')
cmd.color('red', code1 + ' and resi 310')
cmd.color('orange',  code2 + ' and resi 357')
cmd.select('pair_27', '(' + code1 + ' and resi 311) + (' + code2 +' and resi 358)')
cmd.color('red', code1 + ' and resi 311')
cmd.color('orange',  code2 + ' and resi 358')
cmd.select('pair_28', '(' + code1 + ' and resi 312) + (' + code2 +' and resi 359)')
cmd.color('red', code1 + ' and resi 312')
cmd.color('orange',  code2 + ' and resi 359')
cmd.select('pair_29', '(' + code1 + ' and resi 313) + (' + code2 +' and resi 360)')
cmd.color('red', code1 + ' and resi 313')
cmd.color('orange',  code2 + ' and resi 360')
cmd.select('pair_30', '(' + code1 + ' and resi 314) + (' + code2 +' and resi 361)')
cmd.color('red', code1 + ' and resi 314')
cmd.color('orange',  code2 + ' and resi 361')
cmd.select('pair_31', '(' + code1 + ' and resi 315) + (' + code2 +' and resi 362)')
cmd.color('red', code1 + ' and resi 315')
cmd.color('orange',  code2 + ' and resi 362')
cmd.select('pair_32', '(' + code1 + ' and resi 323) + (' + code2 +' and resi 367)')
cmd.color('red', code1 + ' and resi 323')
cmd.color('orange',  code2 + ' and resi 367')
cmd.select('pair_33', '(' + code1 + ' and resi 324) + (' + code2 +' and resi 368)')
cmd.color('red', code1 + ' and resi 324')
cmd.color('orange',  code2 + ' and resi 368')
cmd.select('pair_34', '(' + code1 + ' and resi 325) + (' + code2 +' and resi 369)')
cmd.color('red', code1 + ' and resi 325')
cmd.color('orange',  code2 + ' and resi 369')
cmd.select('pair_35', '(' + code1 + ' and resi 326) + (' + code2 +' and resi 370)')
cmd.color('red', code1 + ' and resi 326')
cmd.color('orange',  code2 + ' and resi 370')
cmd.select('pair_36', '(' + code1 + ' and resi 327) + (' + code2 +' and resi 371)')
cmd.color('red', code1 + ' and resi 327')
cmd.color('orange',  code2 + ' and resi 371')
cmd.select('pair_37', '(' + code1 + ' and resi 328) + (' + code2 +' and resi 372)')
cmd.color('red', code1 + ' and resi 328')
cmd.color('orange',  code2 + ' and resi 372')
cmd.select('pair_38', '(' + code1 + ' and resi 329) + (' + code2 +' and resi 373)')
cmd.color('red', code1 + ' and resi 329')
cmd.color('orange',  code2 + ' and resi 373')
cmd.select('pair_39', '(' + code1 + ' and resi 330) + (' + code2 +' and resi 374)')
cmd.color('red', code1 + ' and resi 330')
cmd.color('orange',  code2 + ' and resi 374')
cmd.select('pair_40', '(' + code1 + ' and resi 331) + (' + code2 +' and resi 375)')
cmd.color('red', code1 + ' and resi 331')
cmd.color('orange',  code2 + ' and resi 375')
cmd.select('pair_41', '(' + code1 + ' and resi 336) + (' + code2 +' and resi 383)')
cmd.color('red', code1 + ' and resi 336')
cmd.color('orange',  code2 + ' and resi 383')
cmd.select('pair_42', '(' + code1 + ' and resi 337) + (' + code2 +' and resi 384)')
cmd.color('red', code1 + ' and resi 337')
cmd.color('orange',  code2 + ' and resi 384')
cmd.select('pair_43', '(' + code1 + ' and resi 338) + (' + code2 +' and resi 385)')
cmd.color('red', code1 + ' and resi 338')
cmd.color('orange',  code2 + ' and resi 385')
cmd.select('pair_44', '(' + code1 + ' and resi 339) + (' + code2 +' and resi 386)')
cmd.color('red', code1 + ' and resi 339')
cmd.color('orange',  code2 + ' and resi 386')
cmd.select('pair_45', '(' + code1 + ' and resi 340) + (' + code2 +' and resi 387)')
cmd.color('red', code1 + ' and resi 340')
cmd.color('orange',  code2 + ' and resi 387')
cmd.select('pair_46', '(' + code1 + ' and resi 341) + (' + code2 +' and resi 388)')
cmd.color('red', code1 + ' and resi 341')
cmd.color('orange',  code2 + ' and resi 388')
cmd.select('pair_47', '(' + code1 + ' and resi 342) + (' + code2 +' and resi 389)')
cmd.color('red', code1 + ' and resi 342')
cmd.color('orange',  code2 + ' and resi 389')
cmd.select('pair_48', '(' + code1 + ' and resi 353) + (' + code2 +' and resi 390)')
cmd.color('red', code1 + ' and resi 353')
cmd.color('orange',  code2 + ' and resi 390')
cmd.select('pair_49', '(' + code1 + ' and resi 354) + (' + code2 +' and resi 391)')
cmd.color('red', code1 + ' and resi 354')
cmd.color('orange',  code2 + ' and resi 391')
cmd.select('pair_50', '(' + code1 + ' and resi 355) + (' + code2 +' and resi 392)')
cmd.color('red', code1 + ' and resi 355')
cmd.color('orange',  code2 + ' and resi 392')
cmd.select('pair_51', '(' + code1 + ' and resi 356) + (' + code2 +' and resi 393)')
cmd.color('red', code1 + ' and resi 356')
cmd.color('orange',  code2 + ' and resi 393')
cmd.select('pair_52', '(' + code1 + ' and resi 357) + (' + code2 +' and resi 394)')
cmd.color('red', code1 + ' and resi 357')
cmd.color('orange',  code2 + ' and resi 394')
cmd.select('pair_53', '(' + code1 + ' and resi 358) + (' + code2 +' and resi 395)')
cmd.color('red', code1 + ' and resi 358')
cmd.color('orange',  code2 + ' and resi 395')
cmd.select('pair_54', '(' + code1 + ' and resi 359) + (' + code2 +' and resi 396)')
cmd.color('red', code1 + ' and resi 359')
cmd.color('orange',  code2 + ' and resi 396')
cmd.select('pair_55', '(' + code1 + ' and resi 360) + (' + code2 +' and resi 397)')
cmd.color('red', code1 + ' and resi 360')
cmd.color('orange',  code2 + ' and resi 397')
cmd.bg_color('gray70')
