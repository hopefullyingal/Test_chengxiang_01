#include "Read.h"
#include "inital.h"
#include "loop.h"


int main() {

	Read mRead;
	Read::data Mdata;
	Mdata = mRead.getValue();

	//mRead.show(Mdata);
	inital mInital;
	inital::inital_run run;
	 
	// �Ȼ�ȡ���ݣ�������
	 mInital.setValue(Mdata);
	 mInital.run();
	 run = mInital.getValue();
	 //mInital.show(run);

	 loop mLoop;
	 loop::newData NEW;

	 mLoop.setValue(run ,Mdata);
	 NEW = mLoop.loop_run();
	 mLoop.draw(NEW);
}