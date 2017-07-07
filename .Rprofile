suppressPackageStartupMessages({
    library(devtools)
    library(testthat)
})

generate_manual <- function(preview=FALSE){
    cmd <- c("R CMD Rd2pdf --force ",
            "--title='FraseR' ",
            "--output='inst/doc/FraseR-manual.pdf'"
    )
    if(!preview){
        cmd <- paste(cmd, "--no-preview", collapse=" ")
    }
    cmd <- paste(cmd, "./", collapse=" ")
    message(cmd)
    system(cmd)
}

install_fraser <- function(){
    document()
    install.packages("./", repos=NULL, type="source")
    generate_manual()
}

load_fraser <- function(){
    document()
    load_all()
}

load_fraser()
