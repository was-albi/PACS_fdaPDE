# INDEX
1. [Intro](#Intro)  
2. [SSH Connectoin](#ssh)   
3. [Basic Commands](#basic_commands) 
	* [Clone a repository](#clone) 
	* [Update a repository](#pull)
	* [Create a repository](#create)
	* [Update a file](#update)
	* [Remove a directory](#remove)
	* [Remove a Repository](#removeRepo)
	* [Rename a directory](#rename)
4. [Branch](#branch)


[git basic](https://git-scm.com/book/it/v1/Basi-di-Git-Repository-Git#Creare-un-repository-in-una-directory-preesistente)

<div id='Intro' \div>

## Intro

### What is git?
Git è diviso in 3 file : 

* directory: working directory
* index: fa da spazio di transito per i files
* HEAD: che punta all'ultimo commit fatto.

Ogni modifica viene VALIDATA con il comando `git commit`, a cui segue il comand `git push` per l'upload nel repository remoto.

Quando creiamo un repository il tutto avviene localmente, lo scopo finale ovviamente è avere un repository remoto.
A tale scopo possiamo sfruttare il social network Github.


### Setta l'identità

```
$ git config --global user.name "John Doe"  
$ git config --global user.email johndoe@example.com
```
controlla le impostazioni:

```
$ git config --list
```
Helper:

```
$ git help <verb>
$ git <verb> --help
$ man git-<verb>
```
[more description here](https://git-scm.com/book/it/v1/Per-Iniziare-Prima-Configurazione-di-Git)

<div id='ssh' \div>

## SSH Connection

##### Create a ssh-key
use `ssh-keygen` command  
For information use `$man ssh-keygen`

Luca Formaggia advice: create rsa and dsa key as well.

```
$ssh-keygen -t dsa -f 'id_dsa'  -C 'myemail@example.com'
$ssh-keygen -t rsa -f 'id_rsa' -C 'myemail@example.com'
```
`-t` specify the type of key, `-C` is the comment, `-f` file name.  
Now these command generate in the ./ssh folder some files:

```
id_dsa	 id_dsa.pub	 id_rsa  id_rsa.pub   know_hosts	
```
The `.pub` files are our pubblic key, which we can share, but the others must be private!

Use `cat id_dsa.pub` to see our key, copy in order to share.

##### Known_hosts file
The `known_hosts`file lets the client authenticate the server, to check that it isn't connecting to an impersonator. The `authorized_keys` file lets the server authenticate the user.

[more on ssh keygen here](https://www.ssh.com/ssh/keygen/)
tele

### GIT connection
[all info](https://help.github.com/articles/connecting-to-github-with-ssh/)

##### 1)Git SSH Connection
1. devi controllare di non avere già una chiave. Di base, le chiavi SSH degli utenti sono salvate nella directory ~/.ssh. Puoi facilmente controllare spostandoti nella directory e controllandone il contenuto:

		$ cd ~/.ssh
		$ ls
		authorized_keys2  id_dsa       known_hosts
		config            id_dsa.pub
		
##### 2)Adding your SSH key to the ssh-agent
[waht is the ssh-agent?](https://www.ssh.com/ssh/agent)

1. Start the ssh-agent in the background.
2. If you're using macOS Sierra 10.12.2 or later, you will need to modify your `~/.ssh/config` file

		Host *
	 	 AddKeysToAgent yes  
	 	 UseKeychain yes
	 	 IdentityFile ~/.ssh/id_rsa 
3. Add your SSH private key to the ssh-agent and store your passphrase in the keychain.
	
		$ ssh-add -K ~/.ssh/id_rsa

NB: Do the same for the dsa key
##### 3)Adding a new SSH key to your GitHub account
1.  Copies the contents of the id_rsa.pub file to your clipboard
		
		$pbcopy < ~/.ssh/id_rsa.pub
2. Go to Github, setting and copy the key

##### 4)Testing your SSH connection
		ssh -T git@github.com

<div id='basic_commands' \div>

## Basic Commands

<div id='clone' \div>

### CLONE a REPOSITORY (via ssh)
Colonare significa "scaricare", piu precisamente è avere una copia in locale di un repository remoto.

Per clonare un repository Git via SSH, puoi specificare un URL ssh:// come questo:

	$ git clone ssh://user@server/project.git
O non specificare proprio il protocollo — Git utilizza SSH non lo specifichi:

	$ git clone user@server:project.git

<div id='pull' \div>

### UPDATE a REPOSITORY
When a remote repository change some file, and we need to update our local reposiroty to last verion (last commit) do a PULL: 

	$git pull 


<div id="create" >

### CREATE a REPOSITORY
Open terminal and go to working directory

1. `$git init` in my folder  

2. Add something to git  
	* `$git add *` for all file and folder in the working directory    
	* `git add myfolder` for a folder 
	* `git add myfile.md` for a single file
3. `$git commit -m "first commit"` per la validazione dell'add   

4. Link our local repository to Github (a remote repository):    
  Go to github and create a new repository "new_repo" WITHOUT README.md

5. Add to remote repository   
`git remote add origin git@github.com:username/new_repo`

6. Update all command in remote repository: push in master branch!  
 `$ git push -u origin master`

to create a new file: `touch README.md`  
[see here simple guide ](http://rogerdudler.github.io/git-guide/index.it.html)

<div id = "update">

### UPDATE a FILE
es: ho modificato il README.md, voglio fare l'upload.  
esegui i tre comandi:

1. `git add README.md` 
2. `git commit -m "upload README.md"` confermo l'add 
3. `git push -u origin master` carico il tutto nel branch master

COmandi utili:

* `git status` per capire a che punto sono del procedimento 
*  `git diff` per vedere quali differenze(modifiche) ci sono tra i file remoti e quelli locali

<div id = "remove">

### REMOVE a DIRECTORY
In the command-line, navigate to your local repository.
Ensure you are in the default branch:
git checkout master

The `$rm -r` command will recursively remove your folder:

1. `$git rm -r folder-name`
2. Commit the change: 
`$git commit -m "Remove duplicated directory"`
3. Push the change to your remote repository:
`$git push origin master`

<div id = "remove">

### REMOVE a DIRECTORY
In the command-line, navigate to your local repository.
Ensure you are in the default branch:
git checkout master

The `$rm -r` command will recursively remove your folder:

1. `$git rm -r folder-name`
2. Commit the change: 
`$git commit -m "Remove duplicated directory"`
3. Push the change to your remote repository:
`$git push origin master`

<div id = "rename">

### RENAMING a REPOSITORY
Open terminal and go to working directory then   
`$git rm -rf .git`

<div id="branch">

## BRANCH

I branch ('ramificazioni') sono utilizzati per sviluppare features che sono isolate l'una dall'altra. Il branch master è quello di default quando crei un repository. Puoi usare altri branch per lo sviluppo ed infine incorporarli ('merge') nel master branch una volta completati.

1. Creo una nuova brach per lavorare sulla beta (codice in working progres)   
`$ git branch beta`  
complimenti ora siamo nel branch 
`beta`  
2. modifichiamo, aggiungiamo ecc nel branch `beta`
3. Per upload in `beta`: To push the current branch and set the remote as upstream, use

 		   git push --set-upstream origin feature_x

4. Volgo unire i 2 repo, passo al master e faccio il merge

		$git checkout master
		$git merge beta
		$git push
		
OSS: per swicciare in vari branch usare:`$git checkout name_branch`  

[More on Branching](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging)
