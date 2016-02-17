proc copybin {newfile dir} {

    if {![file exists $newfile]} {
	error "No such file $newfile"
    }
    if {![file exists $dir]} {
	error "No such directory $dir"
    }
    if {[file isdirectory $newfile]} {
	error "File $newfile is a directory"
    }
    if {![file isdirectory $dir]} {
	error "File $dir is not a directory"
    }

    set filename [file tail $newfile]
    if {[file exists [file join $dir $filename]]} {
	set index 1
	while {[file exists [file join $dir $filename.$index]]} {
	    incr index
	}
	puts "Moving [file join $dir $filename] to [file join $dir $filename.$index]"
	exec mv [file join $dir $filename] [file join $dir $filename.$index]
    }

    puts "Copying $newfile to $dir"
    exec cp $newfile $dir
}

proc copybini {newfile dir} {

    if {![file exists $newfile]} {
	error "No such file $newfile"
    }
    if {![file exists $dir]} {
	error "No such directory $dir"
    }
    if {[file isdirectory $newfile]} {
	error "File $newfile is a directory"
    }
    if {![file isdirectory $dir]} {
	error "File $dir is not a directory"
    }

    exec /usr/sbin/install -f $dir -o $newfile
}

	
	
