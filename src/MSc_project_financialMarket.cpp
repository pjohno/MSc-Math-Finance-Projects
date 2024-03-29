#include "MSc_projects_financialMarket.hpp"
#include <iostream>

namespace MSC_PROJECTS
{
    
    std::ostream& operator<<(std::ostream& output,const Order& O)
    {
        if(O.oType == bid)
            output << " BID ";
        else if(O.oType == ask)
            output << " ASK ";
        output << " t=" << O.day << ":" << O.period;
        output << " p=" << O.price << " v=" << O.vol << " agentID=" << O.agentID;
        return output;
    }
    
    bool operator<(const Order& a,const Order&b){
        if(a.oType < b.oType)return true;
        else if(a.oType==b.oType)
        {
            if(a.oType == bid && a.price>b.price)
            {
                return true;
            }
            else if(a.oType == ask && a.price<b.price)
            {
                return true;
            }
            else if(a.price==b.price)
            {
                if(a.day<b.day)
                    return true;
                else if(a.day==b.day)
                {
                    if(a.period<b.period)
                    {
                        return true;
                    }
                    else 
                        return false;
                }
                else return false;
                
            }
            else
            {
                return false;
            }
        }
        else return false;
    }
    
    std::ostream& operator<<(std::ostream& output,const Agent& A)
    {
        output << " Agent:= " << A.AgentIndex << std::endl;
        output << " Cash:= " << A.cash << " Assets:= "<< A.delta << std::endl;
        output << " Available Cash:= " << A.availableCash << " Available Assets:= "<< A.availableDelta << std::endl;
        //         output << " Order Book:= \n";
        int runningDelta= A.delta;
        int runningCash=A.cash;
        for(const auto& o : A.marketOrders)
        {
            output << " => " << o << "\n";
            if(o.oType == bid)
                runningCash = runningCash - o.vol*o.price;
            else
                runningDelta = runningDelta - o.vol;
        }
        if(fabs(runningCash - A.availableCash)>1.e-8)
            output << "## ERROR :: CASH NOT ADDING UP " << std::endl;
        if(runningDelta!=A.availableDelta)
            output << "## ERROR :: ASSETS NOT ADDING UP " << std::endl;     
        
        output << " Outgoing Link: =>> " << A.outgoingLink << std::endl;
        output << " Incoming Links:";
        for(const auto& l : A.incomingLinks)
        {
            output << " <<= " << l << " ";
        }
        output  << std::endl;
        return output;
    }
    
    void Market::runSimulation(int totalPeriods, int totalIntraPeriods, int tau,double sigma,double A,double sigmaLambda,std::ostream &output)
    {
        // reset the guru at the start of the simulation
        guru=0;
        output << "#period  , #marketPrice  , #guru  , #links, #lambda" << "\n";          
        // run a bunch of 
        for(int period=1;period<=totalPeriods;period++)
        {
            // calculate the current price of the asset from order books
            if(limitOrderBookAsks.size()>0 && limitOrderBookBids.size()>0 )
            {
                double bestBidPrice = limitOrderBookBids.begin()->price;
                double bestAskPrice = limitOrderBookAsks.begin()->price;
                marketPrice = ( bestAskPrice + bestBidPrice ) / 2.;
            }
            
            // create links to another agent
            for(auto& agent : marketAgents)
            {
                // randomly select an agent
                std::uniform_int_distribution<> randomAgent(0,size()-1);
                int possibleOutgoingLink = randomAgent(rng);
                
                // change outgoing link to guru with probability "lambda"
                std::uniform_real_distribution<> probSelectGuru(0.,1.);
                if( probSelectGuru(rng) < lambda )
                    possibleOutgoingLink = guru;
                
                // update links
                if(possibleOutgoingLink!=agent.AgentIndex)
                {
                    // get the agent we used to follow
                    Agent& oldFollowing = marketAgents[agent.outgoingLink];
                    // delete me from their incoming links
                    oldFollowing.incomingLinks.erase(agent.AgentIndex);
                    // now update my link
                    agent.outgoingLink = possibleOutgoingLink;
                    // get the agent we are now following
                    Agent& newFollowing = marketAgents[possibleOutgoingLink];
                    // add me to their incoming links
                    newFollowing.incomingLinks.insert(agent.AgentIndex);
                }
            }
            
            // update agent expectations
            for(auto& agent : marketAgents)
            {
                std::uniform_real_distribution<> sigma0(0.,sigma);
                double percentageLinks = (double)(agent.incomingLinks.size())/(double)(size());
                double networkSigma = sigma0(rng)*(A + percentageLinks*(1.-omega));
                std::normal_distribution<double> eps(0.,networkSigma);
                agent.r = eps(rng);
            }
            
            for(int intraPeriod = 0;intraPeriod<=totalIntraPeriods;intraPeriod++)
            {
                // randomly select an agent
                std::uniform_int_distribution<> U(0,size()-1);
                int randomlySelectedAgent = U(rng);
                Order o = strategy(randomlySelectedAgent,period,intraPeriod);
                submitOrder(o); 
            }

            // find the guru as the agent with the most incoming links
            guru=0;
            unsigned long maxLinks=0;
            for(auto& agent : marketAgents)
            {
                if( agent.incomingLinks.size() > maxLinks )
                {
                    guru = agent.AgentIndex;
                    maxLinks=agent.incomingLinks.size();
                }
            }
            
            // model dlambda / lambda = dW a log normal SDE, bounded by 1 so that lambda <=1
            std::normal_distribution<double> dW(0.,sigmaLambda);
            double phi = dW(rng);
            // if lambda e^phi>1 reflect in the barrier
            if(lambda*exp(phi)>1.)phi=2.*log(1./lambda)-phi;
            lambda = lambda*exp(phi);

            // now delete orders
            deleteOrdersMarket(period-tau);

            output << period << " , " << marketPrice << " , " << guru << " , " << maxLinks << " , " << lambda << "\n";          
            
        }
        output << std::endl;          
              
    }
    
    MSC_PROJECTS::Order Market::strategy(int agentID,int day, int period)
    {
        Order O;
        Agent& trader = marketAgents[agentID];
        // find agent we might follow 
        Agent& following = marketAgents[trader.outgoingLink];
        // calculate price, using trader's "r" and following agent's "r"
        double pt=marketPrice*exp(omega*trader.r+(1-omega)*following.r);
        
        // set time stamp 
        O.day = day;O.period=period;
        // and unique ID
        O.agentID = trader.AgentIndex;
        std::uniform_real_distribution<> U(0.,1.);
        if(pt>marketPrice)
        {
            O.oType = bid;
            O.price = 0;
            double xi1=U(rng),xi2=U(rng);
            int pT = marketPrice*(1-xi1)*(1-xi2) + pt*xi2 + 0.5;
            O.price=std::max(pT,O.price);
            // can only commit to buying given what you've already bud for, so check available cash
            std::uniform_real_distribution<> V(0.,std::max(0.,trader.availableCash)/O.price);
            O.vol = int( V(rng) );
        }
        else
        {
            O.oType = ask;  
            O.price = 2*pt;
            double xi1=U(rng),xi2=U(rng);
            int pT = pt*(1-xi2) + marketPrice*(1+xi1)*xi2 + 0.5;
            O.price=std::min(pT,O.price);
            // can only sell assets you have, so check trader's available assets
            std::uniform_int_distribution<> V(0,std::max(0,trader.availableDelta));
            O.vol = int( V(rng) );
        }
        return O;
        
    }
    
    int Market::submitOrder(MSC_PROJECTS::Order o)
    {
        // get the trader ID from the order
        Agent& trader = marketAgents[o.agentID];
        // check if order already exists, if so exit now
        auto oo = trader.marketOrders.find(o);
        if(oo!=trader.marketOrders.end())return -1;
       
        // first add order to the agents book, 
        if(o.vol>0)
        {
            // adjust available assets at this stage
            if(o.oType== bid)
            {
                double cash=trader.availableCash - o.vol*o.price;
                // if cash not available to execute transaction then exit
                if(cash>=0.)
                    trader.availableCash = cash;
                else
                    return -1;
            }
            else
            {
                double vol= trader.availableDelta - o.vol;
                // if assets not available then exit
                if(vol>=0)
                    trader.availableDelta = vol;
                else 
                    return -1;
            }
            // OK this is valid order
            // put it in the market orders
            trader.marketOrders.insert(o);
        }
        else
        {
            // return from function if agent selects volume = 0  
            return -1;
        }
        
        // next track the order on the market books
        
        // create a copy of the order to track volumes
        Order Onew(o);
        
        // now look at bid offer first
        if(Onew.oType == bid)
        {
            // iterate through the limit order book of ask offers,
            // to see if any match the bid
            // askOIt will be the current best ask 
            auto askOIt = limitOrderBookAsks.begin();
            
            // go through the book
            while( askOIt != limitOrderBookAsks.end() )
            {
                // get a copy of the current best ask order
                Order askCurrentBest(*askOIt);
                // if there are any ask offers below the bid then try and match them
                if( askCurrentBest.price <= Onew.price )
                {
                    // EXECUTE ORDER
                    int volumeRemaining = matchOrder(Onew,askCurrentBest);
                    //                     std::cout << " VOLUME REMAINING " << volumeRemaining << std::endl;
                    if(volumeRemaining<0)
                    {
                        std::cout << "## error matching agents " << std::endl;
                    }
                    
                    // trade executed
                    int volumeExecuted = Onew.vol - volumeRemaining;
                    if(volumeExecuted>0)marketPrice = askCurrentBest.price;
                    // remove executed trades from ask on book
                    askCurrentBest.vol = ( askCurrentBest.vol - volumeExecuted );
                    // remove the ask order from the book
                    askOIt = limitOrderBookAsks.erase(askOIt);
                    
                    // if volume is unexecuted on the ask side,
                    // add it back in
                    if(askCurrentBest.vol>0)
                        limitOrderBookAsks.insert(askCurrentBest);
                    
                    // now check the bid
                    Onew.vol = volumeRemaining;
                    // EXECUTE ORDER FOR BOTH AGENTS
                    if(Onew.vol==0)
                    {
                        // 
                        cancelOrder(o);
                        return 0;
                    }
                }
                else
                {
                    // this mean best ask is above bid, so break from loop
                    break;
                }
            }
            // if we reach this point, it is because there is still some volume left on the bid unmatched
            // if so, add it to the book
            limitOrderBookBids.insert(Onew);
        }
        else
        {
            // iterate through the limit order book of bid offers,
            // to see if any match the ask
            // bidOIt will be the current best ask 
            auto bidOIt = limitOrderBookBids.begin();
            
            // run through bids
            while( bidOIt != limitOrderBookBids.end() )
            {
                Order bidCurrentBest(*bidOIt);
                // if there are any ask offers below the bid then try and match them
                if( bidCurrentBest.price >= Onew.price )
                {
                    // EXECUTE ORDER
                    int volumeRemaining = matchOrder(Onew,bidCurrentBest);
                    if(volumeRemaining<0)
                    {
                        std::cout << "## error matching agents " << std::endl;
                    }
                    
                    // trade executed
                    int volumeExecuted = Onew.vol - volumeRemaining;
                    if(volumeExecuted>0)marketPrice = bidCurrentBest.price;
                    bidCurrentBest.vol = ( bidOIt->vol - volumeExecuted );
                    
                    //                     std::cout << " ERASE OLD BID " << std::endl;
                    bidOIt = limitOrderBookBids.erase(bidOIt);
                    
                    // if there is any volume left on the current best bid, add it back into market
                    if(bidCurrentBest.vol>0)
                    {
                        limitOrderBookBids.insert(bidCurrentBest);
                    }
                    
                    Onew.vol = volumeRemaining;
                    // RETURN 0 if order executed by market book
                    //                     std::cout << " RETURN 0 if order executed by market book " << std::endl;
                    if(Onew.vol==0){
                        cancelOrder(o);
                        return 0;
                    }
                    //                     std::cout << " MOVE TO NEXT BID IN THE BOOK " << std::endl;
                    
                }
                else
                {
                    // 
                    break;
                }
            }
            //             std::cout << " New ask :: " << Onew << std::endl;
            limitOrderBookAsks.insert(Onew);
        }
        return Onew.vol;
    }
    
    int Market::cancelOrder(MSC_PROJECTS::Order O)
    {
        auto& A = marketAgents[O.agentID];
        // check if the order is logged in the agents book
        auto oo = A.marketOrders.find(O);
        // if we can find it remove it
        if(oo != A.marketOrders.end())
        {
            // first update available cash/stock to reflect removal of order
            if(O.oType == bid)
            {
                // order is a bid so we return available cash
                A.availableCash = A.availableCash + O.vol * O.price;
            }
            else
            {
                // order is a ask so we return available stock
                A.availableDelta = A.availableDelta + O.vol;
            }
            // then delete it
            A.marketOrders.erase(O);
            // return 0 if order cancelled by agent
            return 0;
        }
            // return -1 if agent doesn't have this order  
        return -1;
    
    }

    int Market::deleteOrdersMarket(int s)
    {
        auto bidOIt = limitOrderBookBids.begin();
        while(bidOIt != limitOrderBookBids.end())
        {
            if(bidOIt->day < s)
            {
                cancelOrder(*bidOIt);
                bidOIt = limitOrderBookBids.erase(bidOIt);
            }
            else
                bidOIt++;
        }
        auto askOIt = limitOrderBookAsks.begin();
        while(askOIt != limitOrderBookAsks.end())
        {
            if(askOIt->day < s)
            {
                cancelOrder(*askOIt);
                askOIt = limitOrderBookAsks.erase(askOIt);
            }
            else
                askOIt++;
        }
        return 0;
    }

    
    int Market::matchOrder(const MSC_PROJECTS::Order &oExe,const MSC_PROJECTS::Order &oCou)
    {
        // check bid/ask is opposite
        if(oExe.oType == oCou.oType )
        {
            std::cout << "## can't match 2 bids or 2 asks " << std::endl;
            return -1;
        }
        // try to execute v shares at order price
        int unmatchedVol = executeOrder(oCou,oExe.vol,oCou.price);
        if(unmatchedVol<-1)
        {
            std::cout << "## Agent #" << oCou.agentID;
            std::cout << "## has no matching order." << std::endl;
        }
        // volume executed
        int vExecuted = oExe.vol - unmatchedVol;
        // execute executive order
        executeOrder(oExe,vExecuted,oCou.price);
        return unmatchedVol;
    }
    
    
    int Market::executeOrder(const Order &O,int v,int p)
    {
       // get trader
        auto& trader = marketAgents[O.agentID];
        // execute the given order O at the price p, with a volume of at most v
        // check that this agent has the order in his book
        std::set<Order>::iterator oo = trader.marketOrders.find(O);
        // if valid execute order and delete from the agents book
        if(oo != trader.marketOrders.end())
        {
            // try executing minimum volume from Order and requested v
            int volume = std::min(O.vol,v);
//              std::cout << "## EXECUTE " << O << " with vol " << volume << " at price " << p << std::endl;
            // order type
            if(O.oType == bid)
            {
                // check their is sufficient cash to execute. If not execute maximum amount possible.
                volume = std::min(volume,int(trader.cash/p));
                // order is a bid so calculate transaction cost and take it from cash
                double transactionCost = volume * p;
                trader.cash = trader.cash - transactionCost;
                // add stocks 
                trader.delta = trader.delta + volume;
                
                // will be cancelling this order now, so add back in total order cost, and remove actual transaction
                trader.availableCash = trader.availableCash + O.vol * O.price - transactionCost;
                // add stocks to available
                trader.availableDelta = trader.availableDelta + volume;
            }
            else 
            {
                // check delta is sufficiently large to execute. If not execute maximum amount possible.
                volume = std::min(volume,trader.delta);
                // order is an askso calculate transaction cost and take it from cash
                double transactionCost = volume * p;
                trader.cash = trader.cash + transactionCost;
                trader.delta = trader.delta - volume;
                
                // will be cancelling this order now, so add back in total stocks in order, before removing actual amount sold
                trader.availableCash = trader.availableCash + transactionCost;
                trader.availableDelta = trader.availableDelta + O.vol - volume;
            }
            // Create a new order and adjust the volume to the remaining left over after execution
            Order Onew(O);
            Onew.vol = Onew.vol - volume;
            // now remove it from my order list
            trader.marketOrders.erase(oo);
            
            // if any outstanding volume left, add it back in
            // put it in the market orders
            if(Onew.vol>0){
                trader.marketOrders.insert(Onew);
                // adjust available assets at this stage
                if(Onew.oType== bid)
                    trader.availableCash = trader.availableCash - Onew.vol*Onew.price;
                else
                    trader.availableDelta = trader.availableDelta - Onew.vol;
            }
            // return number of orders left unmatched
            return v-volume;
        }
        // return -1 if agent doesn't have this order  
        return -1;
    }
    
    double Agent::wealth(double p) const 
    {
        return cash + delta * p;
    }
    
    void Market::printFullMarketData(std::ostream& output) const
    {
        output << "###############\n# MARKET DATA\n###############\n";
        output << " Market Price " << marketPrice << "\n";
        const Agent& guruAgent = marketAgents[guru];
        output << " Guru " << guruAgent.AgentIndex << " |  Wealth " << guruAgent.wealth(marketPrice) << " | Links " << guruAgent.incomingLinks.size() << "\n";
        output << "#######\n# MARKET BIDS\n######\n";
        for(auto b : limitOrderBookBids)output << b << "\n";
        output << "#######\n# MARKET ASKS\n######\n";
        for(auto a : limitOrderBookAsks)output << a << "\n";
        output << "#######\n# MARKET AGENTS\n######\n";
        for(auto A : marketAgents)output << A << "\n";
    }

}
