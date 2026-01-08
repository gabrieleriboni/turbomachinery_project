# turbomachinery_project

-radiale per uno stadio dato il beta_tt, altrimenti con assiale 2/3 stadi

-radiale perchè industrialmente applicato: non esistono assiali per ammonia, di solito sono volumetrici quindi: dobbiamo scegliere a absso omega_s e quindi quelli con più basssa omega_s sono radiali

-dai calcoli abbiamo un unshrouded poi confermato da calcoli successivi
-da paper china il Set 1 è quello ottimale dato il nostro M1_w_tip
-sistemare solidity con ciclo for, non cambiare thickness
-controllare VISCOSITà DINAMICA IN 1,2,3
-riguardare differenza fra parassitic e internal losses, una va in un efficienza l'altra in un'altra
-trovato problema al codice: viene negativa la T5 perchè troppo alta V5, trovare alternativa