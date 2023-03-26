#!/bin/sh

nm *.o | grep " T " | egrep -v '(parmetis|METIS|__)' | awk '{printf("#define %s libpm4gnn__%s\n",$3, $3)}'

