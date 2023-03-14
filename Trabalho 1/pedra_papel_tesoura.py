from secrets import choice
import numpy as np
import random
print("Vamos jogar pedra, papel e tesoura!!!")

while True:
        print("Primeiro, selecione o modo de jogo entre facil e dificil sem acento")
        dificuldade = input()
        dificuldade = dificuldade.lower()
        if(dificuldade != "dificil" and dificuldade != "facil"):
            print("Selecione entre facil e dificil")
        else:
            break
print("Selecione entre pedra, papel ou tesoura")
while True:
    while True:
        jogador = input()
        jogador = jogador.lower()
        if(jogador != "pedra" and jogador != "papel" and jogador != "tesoura"):
            print("Selecione entre pedra, papel ou tesoura, por favor")
        else:
            break
    lista_aux = [0,1,2]
    if dificuldade == 'dificil':        
        if jogador == 'pedra':
            pesos = [20,20,60]           
        if jogador == 'tesoura':
            pesos = [60,20,20]
        if jogador == 'papel':
            pesos = [20,60,20]
            
        maquina = random.choices(lista_aux, pesos)
        maquina = maquina[0]
    else:
        maquina = np.random.randint(0,3)
   
    if maquina == 0:
        maquina = 'pedra'
    elif maquina == 1:
        maquina = 'tesoura'
    else:
        maquina = 'papel'
        
    if jogador == maquina:
        print('empate. De novo...')
    else:
        if jogador == 'pedra':
            if maquina == 'tesoura':
                print('Você venceu!!!')
            if maquina == 'papel':
                print('Você perdeu. Que pena...')
            
        if jogador == 'papel':
            if maquina == 'pedra':
                print('Você venceu!!!')
            if maquina == 'tesoura':
                print('Você perdeu. Que pena...')
                
        if jogador == 'tesoura':
            if maquina == 'papel':
                print('Você venceu!!!')
            if maquina == 'pedra':
                print('Você perdeu. Que pena...')
        break
    
    