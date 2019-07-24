<!--- Greatly inspired from scikit-learn's contribution guidelines -->

## Contributor's guidelines

This project is a community effort, and everyone is welcome to contribute.

The project is hosted on https://github.com/nservant/HiC-pro


In case you experience issues using this package, do not hesitate to submit a
ticket to the [GitHub issue
tracker](https://github.com/nservant/HiC-Pro/issues). You are also welcome to
post feature requests or pull requests.




### Our community, our values

We are a community based on openness and friendly, didactic, discussions.

We aspire to treat everybody equally, and value their contributions.

Decisions are made based on technical merit and consensus.

Code is not the only way to help the project. Reviewing pull requests,
answering questions to help others on mailing lists or issues, organizing and
teaching tutorials, working on the website, improving the documentation, are
all priceless contributions.

We abide by the principles of openness, respect, and consideration of others
of the [R community](https://user2015.math.aau.dk/behaviouR) 



### Ways to contribute

There are many ways to contribute to HiC-pro, with the most common ones
being contribution of code or documentation to the project. Improving the
documentation is no less important than improving the library itself. If you
find a typo in the documentation, or have made improvements, do not hesitate
to send an email to the mailing list or preferably submit a GitHub pull
request. Full documentation can be found under the doc/ directory.

But there are many other ways to help. In particular answering queries on the
issue tracker, investigating bugs, and reviewing other developers' pull
requests are very valuable contributions that decrease the burden on the
project maintainers.

Another way to contribute is to report issues you're facing, and give a
"thumbs up" on issues that others reported and that are relevant to you. It
also helps us if you spread the word: reference the project from your blog and
articles, link to it from your website, or simply star to say "I use it."

### How to contribute


The preferred way to contribute to HiC-pro is to fork the [main
repository](https://github.com/nservant/HiC-Pro) on GitHub, then submit a
"pull request" (PR).

- [Create an account on GitHub](https://github.com/join) if you do not already
  have one.

- Fork the [project repository](https://github.com/nservant/HiC-pro): click on
  the 'Fork' button near the top of the page. This creates a copy of the code
  under your account on the GitHub user account. For more details on how to
  fork a repository see [this
  guide](https://help.github.com/articles/fork-a-repo/).

- Clone your fork of the HiC-pro repo from your GitHub account to your local disk:

```bash
    $ git clone git@github.com:YourLogin/HiC-pro.git
    $ cd HiC-pro
```

- Create a branch to hold your development changes:

    $ git checkout -b my-feature

    and start making changes. Always use a feature branch. It's good practice
    to never work on the master branch!


- Develop the feature on your feature branch on your computer, using Git to do
  the version control. When you’re done editing, add changed files using git
  add and then git commit files:

```bash
    $ git add modified_files
    $ git commit
```  
    to record your changes in Git, then push the changes to your GitHub account with:

```bash
    $ git push -u origin my-feature
```

    Follow these instructions to create a pull request from your fork. This
    will send an email to the committers. You may want to consider sending an
    email to the mailing list for more visibility.


If any of the above seems like magic to you, then look up the Git
documentation and the Git development workflow on the web, or ask a friend or
another contributor for help.

If some conflicts arise between your branch and the master branch, you need to
merge master. For that, you first need to fetch the upstream, and then merge
its master into your branch:

```bash
$ git fetch upstream
$ git merge upstream/master
```

Subsequently, you need to solve the conflicts. You can refer to the Git
documentation related to resolving merge conflict using the command line.

### Contributing to related projects

HiC-pro thrives in an ecosystem of several related projects, which also may
have relevant issues to work on including:

- [iced](https://github.com/hiclib/iced) 
- [HiTC](https://www.bioconductor.org/packages/release/bioc/html/HiTC.html)
