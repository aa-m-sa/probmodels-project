/*
 *  BEANDisco: interrupt
 *  
 *  Copyright 2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <csignal>
#include <functional>

#ifndef INTERRUPT_HPP
#define INTERRUPT_HPP

/*static struct sigaction act, oldAct;
static bool interruptEnabled = false;

static bool interruptPending = false;

static void signalHandler(int sig) {
	if (sig == SIGINT)
		interruptPending = true;
}

void setInterruptEnabled(bool enable) {
	interruptPending = false;
	if (!interruptEnabled && enable) {
		act.sa_handler = signalHandler;
		sigaction(SIGINT, &act, &oldAct);
		interruptEnabled = true;
	}
	else if (interruptEnabled && !enable) {
		sigaction(SIGINT, &oldAct, &act);
		interruptEnabled = false;
	}
}

bool isInterruptPending() {
	return interruptPending;
}*/

void resetInterruptHandler(std::function<void(int sig)> newHandler = nullptr) {
	static struct sigaction act, oldAct;
	static std::function<void(int sig)> handler = nullptr;

	static auto handlerWrapper = [] (int sig) {
		handler(sig);
	};

	if (newHandler) {
		act.sa_handler = handlerWrapper;
		if (!handler)
			sigaction(SIGINT, &act, &oldAct);
		else
			sigaction(SIGINT, &act, nullptr);
		handler = newHandler;
	}
	else {
		if (handler) {
			sigaction(SIGINT, &oldAct, nullptr);
			handler = nullptr;
		}
	}
}

#endif

