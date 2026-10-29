

### Build vignette
pkgdown::build_articles()

### Build references
pkgdown::build_reference()

### Build home
pkgdown::build_home()


### TO DEPLOY ONLINE
# Preview your site locally before publishing (runs all the previous functions)
pkgdown::build_site()

# Run once to configure your package to use and deploy pkgdown (only works from Baptiste laptop)
usethis::use_github_pages(branch = 'develop', path = "/docs")


