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
	 
	// 先获取数据，再运行
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