package org.broadinstitute.dropseqrna.utils;

public enum Bases { 
	A('A'),
	C('C'),
	G('G'),
	T('T'),
	N('N');
	
	private Character base;
	
	Bases(Character base) {
		this.base=base;
	}
	
	Character getBase() {
		return this.base;
	}
	
}