nome, idade = input('Digite seu nome e idade separados por espaço:').split()
print(f'Seu nome é {nome} e sua idade é {idade} ano(s)')

arquivo = open("amigos.txt", "a")
arquivo.write('Seu nome é ' + nome + ' e sua idade é ' + idade + ' ano(s)\n')