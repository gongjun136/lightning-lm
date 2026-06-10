# 5.5.1.11 fork GitHub 开源框架

## 背景

项目是 GitHub 上的开源框架，但是作者还在不断更新。  
我想基于这个框架添加或修改代码，同时后期作者更新时，我也能比较方便地把最新代码同步过来。

最靠谱的思路其实就一句话：

> 让“上游原版”和“自己的修改版”分开管理。

也就是说，作者的原始仓库作为 `upstream`，我自己的 fork 仓库作为 `origin`，自己的改动放在单独的功能分支里。

---

## 1. 先 fork 原项目

在 GitHub 上点击 **Fork**，得到自己的仓库，例如：

```text
github.com/you/project
```

原作者仓库假设是：

```text
github.com/original/project
```

然后克隆自己的 fork 仓库：

```bash
git clone https://github.com/you/project.git
cd project
```

把原作者仓库添加为 `upstream`：

```bash
git remote add upstream https://github.com/original/project.git
```

查看远程仓库：

```bash
git remote -v
```

一般约定：

```text
origin/main      自己 GitHub 仓库的主分支
upstream/main    原作者仓库的主分支，只用来拉取更新
```

自己的开发不要直接在 `main` 上乱改，尽量新建功能分支。

---

## 2. 自己开发功能

先切到 `main`，并保证它是最新的：

```bash
git checkout main
git pull origin main
```

然后新建自己的功能分支：

```bash
git checkout -b feature/my-awesome-change
```

之后就在这个分支上改代码。

查看修改：

```bash
git status
git diff
```

`git diff` 查看完后，按 `q` 退出。

提交并推送：

```bash
git add .
git commit -m "Add xxx feature"
git push origin feature/my-awesome-change
```

如果自己的改动对原框架也有价值，可以尝试给上游提 PR。  
如果能合进上游，以后这部分就不用自己长期维护了。  
不过现实一点，大佬不一定看得上你的代码。

---

## 3. 上游更新后怎么同步

当原作者仓库更新后，先拉取上游最新代码：

```bash
git fetch upstream
```

然后更新自己仓库的 `main` 分支。

### 方式一：merge

这种方式简单，好理解。

```bash
git checkout main
git merge upstream/main
git push origin main
```

也可以直接写成：

```bash
git checkout main
git pull upstream main
git push origin main
```

### 方式二：rebase

这种方式历史更干净，但对新手来说稍微麻烦一点。

```bash
git checkout main
git rebase upstream/main
git push origin main --force-with-lease
```

如果不熟悉 `rebase`，先用 `merge` 就可以。

---

## 4. 把最新 main 同步到自己的功能分支

比如自己的开发分支是：

```bash
feature/my-awesome-change
```

先切过去：

```bash
git checkout feature/my-awesome-change
```

然后把最新的 `main` 合进来。

### 方式一：merge main

```bash
git merge main
```

### 方式二：rebase main

```bash
git rebase main
```

如果怕麻烦，还是推荐先用 `merge`。  
如果想让提交历史更清爽，再考虑 `rebase`。

---

## 5. 如果出现冲突

如果作者改了你也改过的文件，就可能出现冲突。

可以先把 Git 默认编辑器改成 VS Code，避免进入 Vim：

```bash
git config --global core.editor "code --wait"
```

解决冲突后，把文件重新加入暂存区：

```bash
git add <解决完冲突的文件>
```

如果你是在 `rebase` 过程中出现冲突，继续执行：

```bash
git rebase --continue
```

如果你是在 `merge` 过程中出现冲突，解决后提交即可：

```bash
git commit
```

---

## 简单记法

平时开发：

```bash
git checkout main
git pull origin main
git checkout -b feature/xxx
# 修改代码
git add .
git commit -m "xxx"
git push origin feature/xxx
```

上游更新：

```bash
git fetch upstream
git checkout main
git merge upstream/main
git push origin main
```

把更新同步到自己的开发分支：

```bash
git checkout feature/xxx
git merge main
```

核心原则：

> `main` 尽量保持干净，用来同步上游；自己的改动放到 `feature` 分支里。
