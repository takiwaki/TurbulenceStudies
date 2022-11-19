# GitHubの使い方

## githubとファイルのやりとりができるようにする設定
XCで
    
    ssh-keygen -t rsa
    
ファイル名を聞かれたらid_rsa_gitと答える．  
公開鍵, `id_rsa_git.pub`の中身ををgitに登録する．  

秘密鍵を`./ssh/id_rsa_git`に置く.  

.ssh/configにいかを追加．

    Host github.com
         HostName github.com  
         User     git  
         IdentityFile ~/.ssh/id_rsa_git  
         Port 22
         

## 個人の識別情報
Gitを使いはじめるとき，最初にユーザー名とe-mailアドレスを設定する．
これはマシンを変えるごとに一度設定する必要がある（XCで一回，解析サーバーで一回，京で一回等）．
以下のTaro Tenmonは自分の名前にtaro.tenmonは自分のe-mailアドレスに(GitHubに登録したもの)変更のこと．

    git config --global user.name "Taro Tenmon"
    git config --global user.email taro.tenmon@nao.ac.jp        

## 初回のダウンロード  
`main` ブランチのものを`test1`ディレクトリにダウンロードしたい場合．
    
    git clone -b main git@github.com:takiwaki/TurbulenceStudies test1
    

## アップロード  
アップデートしたら`git status`で変更されたファイルを確認  

`git diff`でどこを変更したのかを確認．  

実際に更新するファイルはaddする（それから先，管理対象になる）  
    
    git add *.F 
    git add Makefile 

全部addしたら  
    
    git commit -m “comment”
    git push origin`  

## 最新版と同期させる
同期させるには以下のコマンドをうつ．
    
    git pull
    
これでうまくマージできない，つまり自分でファイルを変更した分がある場合，それを捨てるのか残すのか決める．まずは`git status`してみる．  
(1)ひとまず退避してとっておく `git stash`  
(2)自分でした変更を削除 `git checkout`  
(3)自分でした変更を残す `git add`   

(1)の場合，`git pull`した後`git stach pop`とする.

## branchをマージする。
    git checkout master
    git merge develope
