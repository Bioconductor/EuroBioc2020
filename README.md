This repository contains material for the European _Bioconductor_
annual conference. [View the conference web site][1].

## Setup

Please follow best practices by previewing changes locally. 

1. Make sure that ruby and bundler are installed, following the
   'Requirements' section of [GitHub's documentation][2]. Use Ruby
   3.6.5.

2. Clone the repository and switch to the `docs/` directory

        git clone git@github.com:Bioconductor/EuroBioc2020
        cd EuroBioc2020/docs

3. Install or update bundler to install the ruby pre-requisities.

        gem install --user-install bundler
        # If the installer complains, add the suggested \$PATH\_TO\_RUBY/bin
        # directory to your ~/.bash_profile or ~/.bashrc or similar.

4. Install ruby pre-requisites.

        bundle config set path 'vendor/bundle'     # once only; references Gemfile
        
5. Update installation according to Gemfile

        bundle install

6. Execute the jekyll server

        bundle exec jekyll serve
        
    and view the results at https://localhost:4000

## Adding content to existing pages

Edit or add material as markdown files in the docs/ directory. Please
wrap lines to 80 character width and aim for simple markdown rather
than elaborate html or other content.

## New pages

Add new markdown files in the `docs/` directory. Files should start with

    ---
    layout: default
    ---
    
    {% include header.md %}

Followed by a level 2 (`##`) heading.

Edit `_include/navigation.html` to link the page to the navigation
sidebar.

[1]: https://bioconductor.github.io/EuroBio2020
[2]: https://help.github.com/articles/setting-up-your-github-pages-site-locally-with-jekyll/#requirements
