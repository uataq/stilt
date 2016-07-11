#---------------------------------------------------------------------------------------------------
#  $Id: unix.shell.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

unix.shell<-function(command, shell = "/bin/sh", ...){
        tfile <- tempfile("Sshell")
        on.exit(unlink(tfile))
        cat(file = tfile, command, "\n", sep = "")
        command <- paste(shell, tfile)
        unix(command, ...)
}
