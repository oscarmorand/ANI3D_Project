#pragma once

namespace cgp
{

	class timer_basic
	{
	public:
		timer_basic();
		float update();
		void start();
		void stop();
		bool is_running() const;

		float t;
		float scale;

	protected:
		bool running;
		float time_previous;
	};

}