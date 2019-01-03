#! /usr/bin/env python
# -*- coding:utf-8 -*-
import smtplib
from email.mime.text import MIMEText
from email.utils import formataddr
import smtplib
from email.mime.text import MIMEText
from email.utils import formataddr


def sendemail():
    try:
        msg=MIMEText('Notification!Done!','plain','utf-8')
        msg['From']=formataddr(["R",'chency1997@126.com'])
        msg['To']=formataddr(["ccy",'710969718@qq.com'])
        msg['Subject']="Done!"
        
        server=smtplib.SMTP_SSL("smtp.126.com",465)
        server.login("chency1997@126.com","2010ccy")
        server.sendmail('chency1997@126.com',['710969718@qq.com','cychenbnu@icloud.com','chency1997@126.com'],msg.as_string())
        server.quit()
        print "Succeed"
    except:
        print "False"

sendemail()
